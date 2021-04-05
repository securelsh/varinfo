#include "../header/varinfo.h"
#include "../header/statisticaltest.h"
#include <cmath>


pthread_mutex_t mutex_filter = PTHREAD_MUTEX_INITIALIZER;

using namespace std;


/*

- To Do:
  -- change variable name front strand to forward strand.
  
- Changelog:
  -- Date Mar-31-2021 SeokCholHong shulkhorn@gmail.com
     --- changed ref/alt conditions in contingency table. alt is counted based on adi alt

  -- Date Mar-24-2021 SeokCholHong shulkhorn@gmail.com
     --- fill variant info to VARIANT struct

  -- Date Mar-22-2021 SeokCholHong shulkhorn@gmail.com
     --- make method for estimation of strand bias & VAF

  -- Date Mar-16-2021 SeokCholHong shulkhorn@gmail.com
     --- Extract information from BAM file

*/


bool CINFO::CalcInfo()
{
	if(m_bIsDebug) cout << "- CINFO::CalcInfo()" << endl;

	m_Input.vfVaf.resize(m_Input.vsChr.size());
	m_Input.vfStrandBias.resize(m_Input.vsChr.size());
	m_Input.v2nReadLen.resize(m_Input.vsChr.size());
	m_Input.v2nMapQ.resize(m_Input.vsChr.size());
	m_Input.v2nBaseQ.resize(m_Input.vsChr.size());
	
	vector<DIST_THREAD> vDistThread;
	for(int i=0; i<m_nCntThread; i++)
	{   
		DIST_THREAD DistThread;
		DistThread.INFO = this;
		DistThread.nId = i;
		vDistThread.push_back(DistThread);
	}

	this->SetStartTime();

	pthread_t* thread_handles;
	thread_handles = new pthread_t[m_nCntThread];
	for(int i=0; i<m_nCntThread; i++)	pthread_create(&thread_handles[i], NULL, AlleleDist_helper, (void*)&vDistThread[i]);
	for(int i=0; i<m_nCntThread; i++)	pthread_join(thread_handles[i], NULL);
	free(thread_handles);
	cout << endl;

	return true;
}


void* CINFO::AlleleDist(int nId)
{
	bamFile finBam;
	finBam = bam_open(m_sBamFile.c_str(), "r");
	if(finBam == NULL)	throw std::logic_error("ERROR: Can not open .bam file");

	bam_header_t *bamHeader;
	bamHeader = bam_header_read(finBam);
	if(!bamHeader)       throw std::logic_error("ERROR: Fail to read the header of a bam file");

	bam_index_t *bamIndex;
	bamIndex = bam_index_load(m_sBamFile.c_str());
	if(!bamIndex)      throw std::logic_error("ERROR: Fail to read .bai file");

	int nSubCnt = (int)m_Input.vsChr.size()/m_nCntThread;
	int nExtra = (int)m_Input.vsChr.size()%m_nCntThread;
	if(nId < nExtra)	nSubCnt++;

	for(int i=nId; i<(int)m_Input.vsChr.size(); i+=m_nCntThread)
	{
		string sChr, sRef, sAlt;
		int nChr, nPos;
		
		sChr = m_Input.vsChr[i];
		sRef = m_Input.vsRef[i];
		sAlt = m_Input.vsAlt[i];
		nPos = m_Input.vnPos[i];
		nChr = ConvertChrToTid(sChr, bamHeader);
		if(nChr == -1)
		{
			cout << sChr << "(" << nChr << ")" << ":" << nPos << "\t" << sRef << "\t" << sAlt << endl;
			throw std::logic_error("Not available Chr");
		}
		
		this->GetVarInfo(i, sChr, nChr, nPos, sRef, sAlt, bamIndex, finBam);

	}

	bam_close(finBam);
	
	return (void*)1;

}

bool CINFO::GetVarInfo(int nIdx, string sChr, int nChr, int nPos, string sRef, string sAlt, bam_index_t *bamIndex, bamFile finBam)
{
	bam_iter_t bamIter;
	bam1_t *b;
	b = bam_init1();
	bamIter = bam_iter_query(bamIndex, nChr, nPos-1, nPos);

	QCINFO QcInfo;
	int nRet;
	while((nRet = bam_iter_read(finBam, bamIter, b)) >= 0){
		GetNb(b, nPos - (b->core.pos+1), QcInfo);
	}
	GetAnalysis(QcInfo, sRef, sAlt);

	m_Input.vfVaf[nIdx] = QcInfo.fVaf;
	m_Input.vfStrandBias[nIdx] = QcInfo.fStrandBias;
	m_Input.v2nReadLen[nIdx] = QcInfo.vnReadLen;
	m_Input.v2nMapQ[nIdx] = QcInfo.vnMapQ;
	m_Input.v2nBaseQ[nIdx] = QcInfo.vnBaseQ;

	if(QcInfo.vnMapQ.size()<1||QcInfo.vnBaseQ.size()<1){
		cout << sChr << " " << nPos << endl;
		exit(1);
	}
	if(m_Input.v2nMapQ[nIdx].size()<1||m_Input.v2nBaseQ[nIdx].size()<1){
		cout << sChr << " " << nPos << endl;
		exit(1);
	}

	bam_iter_destroy(bamIter);
	bam_destroy1(b);

	return false;
}

bool CINFO::GetNb(bam1_t *b, int nPos, QCINFO &BamInfo)
{
	bool bIsFit = false;
	uint32_t nTemp=0;
	int l,nPl,nPreAlgn;
	for(l=b->core.l_qname, nPl=1, nPreAlgn = 0;
		l<(b->core.l_qname + (b->core.n_cigar*4));
		l++, nPl++){
		int nChk = nPl%4;
		if(nChk==0){
			nTemp = (b->data[l]<<24) | nTemp;

			int nAlgn = (nTemp>>4);
			int op = nTemp&0xf;
			if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
				nPreAlgn += nAlgn;
				if(nPos<nPreAlgn){
					if(nPos==nPreAlgn-1) bIsFit = true;
					break;
				}
			}
			else if(op == BAM_CINS){
				if(nPos<nPreAlgn){
					if(nPos==nPreAlgn-1) bIsFit = true;
					break;
				}
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CDEL){
				if(nPreAlgn<=nPos && nPos<nPreAlgn+nAlgn){
					return true;
				}
				else if(nPos>=nPreAlgn+nAlgn){
					nPos -= nAlgn;
				}
				else if(nPos<nPreAlgn){
					if(nPos==nPreAlgn-1) bIsFit = true;
					break;
				}
			}
			else if(op == BAM_CSOFT_CLIP || op == BAM_CPAD){
				if(nPos<nPreAlgn){
					if(nPos==nPreAlgn-1) bIsFit = true;
					break;
				}
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CREF_SKIP){
				if(nPos<nPreAlgn){
					if(nPos==nPreAlgn-1) bIsFit = true;
					break;
				}
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CHARD_CLIP){
			}
		}
		else if(nChk==1)	nTemp = b->data[l];
		else if(nChk==2)	nTemp = (b->data[l]<<8) | nTemp;
		else if(nChk==3)	nTemp = (b->data[l]<<16) | nTemp;
	}
	if(nPos>=b->core.l_qseq)	return 0;

	bool bIsIns=false;
	bool bIsDel=false;
	int nIndelLen;
	if(bIsFit){		
		while(l<(b->core.l_qname + (b->core.n_cigar*4))){
			l++, nPl++;
			int nChk = nPl%4;
			if(nChk==0){
				nTemp = (b->data[l]<<24) | nTemp;
				nIndelLen = (nTemp>>4);
				int op = nTemp&0xf;
				if(op == BAM_CINS){
					bIsIns=true;
				}
				else if(op == BAM_CDEL){
					bIsDel=true;
				}
				break;
			}
			else if(nChk==1)	nTemp = b->data[l];
			else if(nChk==2)	nTemp = (b->data[l]<<8) | nTemp;
			else if(nChk==3)	nTemp = (b->data[l]<<16) | nTemp;
		}
	}

	if(bam1_strand(b))
	{
		if(bIsIns){
			string sInsSeq="";
			for(int len=0;len<=nIndelLen;len++){
				unsigned char ucBase = bam1_seqi(bam1_seq(b), nPos+len);
				if     (ucBase == 0x01)	sInsSeq+="A";
				else if(ucBase == 0x02)	sInsSeq+="C";
				else if(ucBase == 0x04)	sInsSeq+="G";
				else if(ucBase == 0x08)	sInsSeq+="T";
				else if(ucBase == 0x0F)	sInsSeq+="N";
				else
					throw std::logic_error("ERROR: sequence error. It deosn't belong to ACGTN (in function CCOMPDP::GetNb():64)");
			}
			BamInfo.vsReverseIns.push_back(sInsSeq);
		} else if(bIsDel) {
			BamInfo.vnReverseDel.push_back(nIndelLen);
		} else {
			unsigned char ucBase = bam1_seqi(bam1_seq(b), nPos);
			if     (ucBase == 0x01)	BamInfo.nReverseA++;
			else if(ucBase == 0x02)	BamInfo.nReverseC++;
			else if(ucBase == 0x04)	BamInfo.nReverseG++;
			else if(ucBase == 0x08)	BamInfo.nReverseT++;
			else if(ucBase == 0x0F)	return 0;
			else
				throw std::logic_error("ERROR: sequence error. It deosn't belong to ACGTN (in function CCOMPDP::GetNb():64)");
		}
	}
	else
	{
		if(bIsIns){
			string sInsSeq="";
			for(int len=0;len<=nIndelLen;len++){
				unsigned char ucBase = bam1_seqi(bam1_seq(b), nPos+len);
				if     (ucBase == 0x01)	sInsSeq+="A";
				else if(ucBase == 0x02)	sInsSeq+="C";
				else if(ucBase == 0x04)	sInsSeq+="G";
				else if(ucBase == 0x08)	sInsSeq+="T";
				else if(ucBase == 0x0F)	sInsSeq+="N";
				else
					throw std::logic_error("ERROR: sequence error. It deosn't belong to ACGTN (in function CCOMPDP::GetNb():64)");
			}
			BamInfo.vsForwardIns.push_back(sInsSeq);
		} else if(bIsDel) {
			BamInfo.vnForwardDel.push_back(nIndelLen);
		} else {
			unsigned char ucBase = bam1_seqi(bam1_seq(b), nPos);
			if     (ucBase == 0x01)	BamInfo.nFrontA++;
			else if(ucBase == 0x02)	BamInfo.nFrontC++;
			else if(ucBase == 0x04)	BamInfo.nFrontG++;
			else if(ucBase == 0x08)	BamInfo.nFrontT++;
			else if(ucBase == 0x0F)	return 0;
			else
				throw std::logic_error("ERROR: sequence error. It deosn't belong to ACGTN (in function CCOMPDP::GetNb():64)");
		}
	}
	
	BamInfo.vnReadLen.push_back(b->core.l_qseq);
	char *nBaseQual = (char*)calloc(b->core.l_qseq, 1);
	memcpy(nBaseQual, bam1_qual(b), b->core.l_qseq);
	BamInfo.vnBaseQ.push_back(nBaseQual[nPos]);
	BamInfo.vnMapQ.push_back(b->core.qual);

	return 1;
}

bool CINFO::GetAnalysis(QCINFO &QcInfo, string sRef, string sAlt){
	int nFR = 0;
	int nFA = 0;
	int nRR = 0;
	int nRA = 0;

	uint16_t nFrontA       = QcInfo.nFrontA    ;
	uint16_t nFrontC       = QcInfo.nFrontC    ;
	uint16_t nFrontG       = QcInfo.nFrontG    ;
	uint16_t nFrontT       = QcInfo.nFrontT    ;
	uint16_t nForwardIns   = QcInfo.vsForwardIns.size();
	uint16_t nForwardDel   = QcInfo.vnForwardDel.size();
	uint16_t nReverseA     = QcInfo.nReverseA  ;
	uint16_t nReverseC     = QcInfo.nReverseC  ;
	uint16_t nReverseG     = QcInfo.nReverseG  ;
	uint16_t nReverseT     = QcInfo.nReverseT  ;
	uint16_t nReverseIns   = QcInfo.vsReverseIns.size();
	uint16_t nReverseDel   = QcInfo.vnReverseDel.size();

	int nRefSize=sRef.size();
	int nAltSize=sAlt.size();
	if(nRefSize<nAltSize){
		int nMatch=0;
		for(int i=0;i<nForwardIns;i++){
			if(QcInfo.vsForwardIns[i]==sAlt)	nMatch++;
		}
		nFR = nFrontA+nFrontC+nFrontG+nFrontT+nForwardDel+(nForwardIns-nMatch);
		nFA = nMatch;
		nMatch=0;
		for(int i=0;i<nReverseIns;i++){
			if(QcInfo.vsReverseIns[i]==sAlt)	nMatch++;
		}
		nRR = nReverseA+nReverseC+nReverseG+nReverseT+nReverseDel+(nReverseIns-nMatch);
		nRA = nMatch;

		QcInfo.fStrandBias = STATTEST::GetFisherPvalue(nFR,nRR,nFA,nRA);
	}
	else if(nRefSize>nAltSize){
		int nMatch=0;
		for(int i=0;i<nForwardDel;i++){
			if(QcInfo.vnForwardDel[i]==(nRefSize-nAltSize))		nMatch++;
		}
		nFR = nFrontA+nFrontC+nFrontG+nFrontT+nForwardIns+(nForwardDel-nMatch);
		nFA = nMatch;
		nMatch=0;
		for(int i=0;i<nReverseDel;i++){
			if(QcInfo.vnReverseDel[i]==(nRefSize-nAltSize))		nMatch++;
		}
		nRR = nReverseA+nReverseC+nReverseG+nReverseT+nReverseIns+(nReverseDel-nMatch);
		nRA = nMatch;

		QcInfo.fStrandBias = STATTEST::GetFisherPvalue(nFR,nRR,nFA,nRA);
	}
	else{
		if(sAlt=="A"){
			nFA = nFrontA;
			nFR = nFrontC+nFrontG+nFrontT+nForwardIns+nForwardDel;
			nRA = nReverseA;
			nRR = nReverseC+nReverseG+nReverseT+nReverseIns+nReverseDel;
		}
		else if(sAlt=="C"){
			nFA = nFrontC;
			nFR = nFrontA+nFrontG+nFrontT+nForwardIns+nForwardDel;
			nRA = nReverseC;
			nRR = nReverseA+nReverseG+nReverseT+nReverseIns+nReverseDel;
		}
		else if(sAlt=="G"){
			nFA = nFrontG;
			nFR = nFrontC+nFrontA+nFrontT+nForwardIns+nForwardDel;
			nRA = nReverseG;
			nRR = nReverseC+nReverseA+nReverseT+nReverseIns+nReverseDel;
		}
		else if(sAlt=="T"){
			nFA = nFrontT;
			nFR = nFrontC+nFrontG+nFrontA+nForwardIns+nForwardDel;
			nRA = nReverseT;
			nRR = nReverseC+nReverseG+nReverseA+nReverseIns+nReverseDel;
		}
		QcInfo.fStrandBias = STATTEST::GetFisherPvalue(nFR,nRR,nFA,nRA);
		QcInfo.fVaf = (float)(nFA+nRA)/(nFR+nRR+nFA+nRA);
	}

	if(nRefSize!=nAltSize){
		if(sRef[0]=='A'){
			nFR = nFrontA;
			nFA = nFrontC+nFrontG+nFrontT+nForwardIns+nForwardDel;
			nRR = nReverseA;
			nRA = nReverseC+nReverseG+nReverseT+nReverseIns+nReverseDel;
		}
		else if(sRef[0]=='C'){
			nFR = nFrontC;
			nFA = nFrontA+nFrontG+nFrontT+nForwardIns+nForwardDel;
			nRR = nReverseC;
			nRA = nReverseA+nReverseG+nReverseT+nReverseIns+nReverseDel;
		}
		else if(sRef[0]=='G'){
			nFR = nFrontG;
			nFA = nFrontC+nFrontA+nFrontT+nForwardIns+nForwardDel;
			nRR = nReverseG;
			nRA = nReverseC+nReverseA+nReverseT+nReverseIns+nReverseDel;
		}
		else if(sRef[0]=='T'){
			nFR = nFrontT;
			nFA = nFrontC+nFrontG+nFrontA+nForwardIns+nForwardDel;
			nRR = nReverseT;
			nRA = nReverseC+nReverseG+nReverseA+nReverseIns+nReverseDel;
		}
		QcInfo.fVaf = (float)(nFA+nRA)/(nFR+nRR+nFA+nRA);
	}
	return true;
}

int CINFO::ConvertChrToTid(string sChr, bam_header_t* bamHeader)
{
	int nId = -1;
	for(int i=0; i<(int)bamHeader->n_targets; i++)
	{
		if(strcmp(sChr.c_str(), bamHeader->target_name[i]) == 0)		
		{
			nId=i;
			break;
		}
	}
	return nId;
}
















