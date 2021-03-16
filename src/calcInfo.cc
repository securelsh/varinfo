#include "../header/varinfo.h"
#include <math.h>


pthread_mutex_t mutex_filter = PTHREAD_MUTEX_INITIALIZER;

using namespace std;


/*

- To Do:
  -- make method for estimation of strand bias & VAF
  -- fill variant info to VARIANT struct
  
- Changelog:
  -- Date Mar-16-2021 SeokCholHong shulkhorn@gmail.com
     --- Extract information from BAM file

*/



//get information
bool CINFO::CalcInfo()
{
	if(m_bIsDebug) cout << "- CINFO::CalcInfo()" << endl;
	
	//parallel process
	vector<DIST_THREAD> vDistThread;
	for(int i=0; i<m_nCntThread; i++)
	{   
		DIST_THREAD DistThread;
		DistThread.INFO = this;
		DistThread.nId = i;
		vDistThread.push_back(DistThread);
	}

	this->SetStartTime();

	//create thread
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
	//open bam
	bamFile finBam;
	finBam = bam_open(m_sBamFile.c_str(), "r");
	if(finBam == NULL)	throw std::logic_error("ERROR: Can not open .bam file");

	bam_header_t *bamHeader;
	bamHeader = bam_header_read(finBam);
	if(!bamHeader)       throw std::logic_error("ERROR: Fail to read the header of a bam file");

	bam_index_t *bamIndex;
	bamIndex = bam_index_load(m_sBamFile.c_str());
	if(!bamIndex)      throw std::logic_error("ERROR: Fail to read .bai file");

	//calc info
	int nSubCnt = (int)m_Input.vsChr.size()/m_nCntThread;
	int nExtra = (int)m_Input.vsChr.size()%m_nCntThread;
	if(nId < nExtra)	nSubCnt++;

	//TODO: allocate a structure
	


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
		
		//Main function
		this->GetVarInfo(i, sChr, nChr, nPos, sRef, sAlt, bamIndex, finBam);

	}

	bam_close(finBam);
	
	return (void*)1;

}

bool CINFO::GetVarInfo(int nIdx, string sChr, int nChr, int nPos, string sRef, string sAlt, bam_index_t *bamIndex, bamFile finBam)
{
	// read bam file
	bam_iter_t bamIter;
	bam1_t *b;
	b = bam_init1();
	bamIter = bam_iter_query(bamIndex, nChr, nPos-1, nPos);// nSPos<= locus <=nEPos

	QCINFO QcInfo;
	int nRet;
	while((nRet = bam_iter_read(finBam, bamIter, b)) >= 0){
		GetNb(b, nPos - (b->core.pos+1), QcInfo);
	}
	GetStrandBias(QcInfo);
	GetVaf(QcInfo);

	// if(QcInfo.nDepthIndel>0){
	// 	cout << sChr << " " << nPos << " " << sRef << " " << sAlt << endl;
	// 	cout << QcInfo.nReverse << " " << QcInfo.nDepthIndel << " " << QcInfo.nDepthA << " " << QcInfo.nDepthC << " " << QcInfo.nDepthG << " " << QcInfo.nDepthT << endl;
	// 	for(unsigned int q=0;q<QcInfo.nMapQ.size();q++)
	// 		cout << (int)QcInfo.nMapQ[q] << " ";
	// 	cout << endl;
	// 	for(unsigned int q=0;q<QcInfo.nBaseQ.size();q++)
	// 		cout << (char)(QcInfo.nBaseQ[q]+33);
	// 	cout << endl;
	// 	exit(1);
	// }

	bam_iter_destroy(bamIter);
	bam_destroy1(b);

	return false;
}

bool CINFO::GetNb(bam1_t *b, int nPos, QCINFO &BamInfo)
{
	//cigar
	uint32_t nTemp=0;
	int l,nPl,nPreAlgn;
	for(l=b->core.l_qname, nPl=1, nPreAlgn = 0;
		l<(b->core.l_qname + (b->core.n_cigar*4));
		l++, nPl++){
		int nChk = nPl%4;
		if(nChk==0){
			nTemp = (b->data[l]<<24) | nTemp;

			int nAlgn = (nTemp>>4);//the higher 28 bits keep the length of a CIGAR
			int op = nTemp&0xf;//the lower 4 bits gives a CIGAR operation
			if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
				nPreAlgn += nAlgn;
				if(nPos<nPreAlgn)	break;
			}
			else if(op == BAM_CINS){//insertion
				if(nPos<nPreAlgn)	break;
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CDEL){//deletion
				if(nPreAlgn<=nPos && nPos<nPreAlgn+nAlgn){
					return true;
				}
				else if(nPos>=nPreAlgn+nAlgn){
					nPos -= nAlgn;
				}
				else if(nPos<nPreAlgn)	break;
			}
			else if(op == BAM_CSOFT_CLIP || op == BAM_CPAD){
				//soft clip|padding
				if(nPos<nPreAlgn)	break;
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CREF_SKIP){//skipping|spliced
				if(nPos<nPreAlgn)	break;
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CHARD_CLIP){
				//hard clipping
			}
		}
		else if(nChk==1)	nTemp = b->data[l];
		else if(nChk==2)	nTemp = (b->data[l]<<8) | nTemp;
		else if(nChk==3)	nTemp = (b->data[l]<<16) | nTemp;
	}
	if(nPos>=b->core.l_qseq)	return 0;

	// check indel
	bool bIsIndel=false;
	while(l<(b->core.l_qname + (b->core.n_cigar*4))){
		l++, nPl++;
		int nChk = nPl%4;
		if(nChk==0){
			nTemp = (b->data[l]<<24) | nTemp;
			int op = nTemp&0xf;//the lower 4 bits gives a CIGAR operation
			if(op == BAM_CINS || op == BAM_CDEL){
				bIsIndel=true;
			}
		}
		else if(nChk==1)	nTemp = b->data[l];
		else if(nChk==2)	nTemp = (b->data[l]<<8) | nTemp;
		else if(nChk==3)	nTemp = (b->data[l]<<16) | nTemp;
	}

	if(bIsIndel) {
		BamInfo.nDepthIndel++;
	} else {
		// unsigned char ucBase = (b->data[(b->core.l_qname + (b->core.n_cigar*4))+(nPos>>1)] >> ((~nPos&1)<<2) & 0xf);
		unsigned char ucBase = bam1_seqi(bam1_seq(b), nPos);
		if     (ucBase == 0x01)	BamInfo.nDepthA++;
		else if(ucBase == 0x02)	BamInfo.nDepthC++;
		else if(ucBase == 0x04)	BamInfo.nDepthG++;
		else if(ucBase == 0x08)	BamInfo.nDepthT++;
		else if(ucBase == 0x0F)	return 0; //N : missing, no depth
		else
			throw std::logic_error("ERROR: sequence error. It deosn't belong to ACGTN (in function CCOMPDP::GetNb():64)");
	}
	// base quality
	char *nBaseQual = (char*)calloc(b->core.l_qseq, 1);
	memcpy(nBaseQual, bam1_qual(b), b->core.l_qseq);
	BamInfo.nBaseQ.push_back(nBaseQual[nPos]);
	// for(unsigned int qq=0;qq<b->core.l_qseq;qq++)
	// 	cout << (char)(nBaseQual[qq]+33);
	// cout << endl;

	// mapping quality
	BamInfo.nMapQ.push_back(b->core.qual);

	// number of reverse strand
	BamInfo.nReverse+=bam1_strand(b);

	// // depth
	// BamInfo.nDepth+=1;

	return 1;
}

/*
 -- Reference
	Guo, Y., Li, J., Li, CI. et al. 
	The effect of strand bias in Illumina short-read sequencing data.
	BMC Genomics 13, 666 (2012). 
	https://doi.org/10.1186/1471-2164-13-666

	< Table 1 >
	https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-666/tables/1
*/
bool CINFO::GetStrandBias(QCINFO &QcInfo){
	float fStrandBias;
	int nFR = 0; // Forward strand reference allele.
	int nFA = 0; // Forward strand non reference allele.
	int nRR = 0; // Reverse strand reference allele.
	int nRA = 0; // Reverse strand non reference allele.

	QcInfo.fStrandBias = fStrandBias;

	return true;
}

bool CINFO::GetVaf(QCINFO &QcInfo){
	float fVaf;

	QcInfo.fVaf = fVaf;

	return true;
}

/*

bool CRD::FilterRead(string sChr, vector<int> &vnNorSPos, vector<int> &vnNorEPos, 
		vector<int> &vnVarSPos, vector<int> &vnVarEPos,
		vector<string> &vsNorSeq, vector<string> &vsVarSeq, int &nSumNor, int &nSumVar, int &nCntNor, int &nCntVar)
{
	vector<int> vnNorSPosT = vnNorSPos;
	vector<int> vnNorEPosT = vnNorEPos;
	vector<int> vnVarSPosT = vnVarSPos;
	vector<int> vnVarEPosT = vnVarEPos;
	vector<string> vsNorSeqT = vsNorSeq;
	vector<string> vsVarSeqT = vsVarSeq;
	vnNorSPos.clear();
	vnNorEPos.clear();
	vnVarSPos.clear();
	vnVarEPos.clear();
	vsNorSeq.clear();
	vsVarSeq.clear();
	
	//calculate midium read length
	vector<int> vnReadLen;
	int nMidLen = 0;
	for(int i=0; i<vnNorSPosT.size(); i++)		vnReadLen.push_back(vnNorEPosT[i]-vnNorSPosT[i]+1);
	for(int i=0; i<vnVarSPosT.size(); i++)		vnReadLen.push_back(vnVarEPosT[i]-vnVarSPosT[i]+1);
	sort(vnReadLen.begin(), vnReadLen.end());
	if(vnReadLen.size() != 0)					nMidLen = vnReadLen[vnReadLen.size()/2];
	
	//calculate # of mutations
	int nNorMidMutation = 0;
	int nVarMidMutation = 0;
	int nPivotSPos = -1;
	int nPivotEPos = -1;
	for(int i=0; i<vnNorSPosT.size(); i++)
	{
		if(nPivotSPos == -1 || vnNorSPos[i] < nPivotSPos)	nPivotSPos = vnNorSPos[i];
		if(nPivotEPos == -1 || vnNorEPos[i] > nPivotEPos)	nPivotEPos = vnNorEPos[i];
	}
	for(int i=0; i<vnVarSPosT.size(); i++)
	{
		if(nPivotSPos == -1 || vnVarSPos[i] < nPivotSPos)	nPivotSPos = vnVarSPos[i];
		if(nPivotEPos == -1 || vnVarEPos[i] > nPivotEPos)	nPivotEPos = vnVarEPos[i];
	}

	string sRefSeq;
	if(!m_FaFile.GetSeq(sRefSeq, sChr, nPivotSPos, nPivotEPos))		return false;
	
	vector<int> vnNorMismatch;
	vector<int> vnVarMismatch;
	
	for(int i=0; i<vnNorSPosT.size(); i++)
	{
		int nCntMismatch = 0;
		for(int j=vnNorSPosT[i]; j<=vnNorEPosT[i]; j++)
			if(vsNorSeqT[i][j-vnNorSPosT[i]] != sRefSeq[j-nPivotSPos])		nCntMismatch++;
		vnNorMismatch.push_back(nCntMismatch);
	}
	for(int i=0; i<vnVarSPosT.size(); i++)
	{
		int nCntMismatch = 0;
		for(int j=vnVarSPosT[i]; j<=vnVarEPosT[i]; j++)
			if(vsVarSeqT[i][j-vnVarSPosT[i]] != sRefSeq[j-nPivotSPos])		nCntMismatch++;
		vnVarMismatch.push_back(nCntMismatch);
	}
	int nMidMutation = 5;
	sort(vnVarMismatch.begin(), vnVarMismatch.end());
	sort(vnNorMismatch.begin(), vnNorMismatch.end());
	if(vnNorMismatch.size() > vnVarMismatch.size())		nMidMutation = vnNorMismatch[vnNorMismatch.size()/2];
	else												nMidMutation = vnVarMismatch[vnVarMismatch.size()/2];

	
	
	// read filtering
	int pnVarA[500] = {0,};
	int pnVarC[500] = {0,}; 
	int pnVarG[500] = {0,}; 
	int pnVarT[500] = {0,}; 
	int pnVarAll[500] = {0,};
	int pnNorA[500] = {0,};
	int pnNorC[500] = {0,}; 
	int pnNorG[500] = {0,}; 
	int pnNorT[500] = {0,}; 
	int pnNorAll[500] = {0,};


	for(int i=0; i<vnVarSPosT.size(); i++)
	{
		for(int j=vnVarSPosT[i]; j<=vnVarEPosT[i]; j++)
		{
			int nPos = j-nPivotSPos;
			if(nPos > 500)	continue;
			if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'A')		pnVarA[nPos]++;
			if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'C')		pnVarC[nPos]++;
			if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'G')		pnVarG[nPos]++;
			if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'T')		pnVarT[nPos]++;
			pnVarAll[nPos]++;
		}
	}
	for(int i=0; i<vnNorSPosT.size(); i++)
	{
		for(int j=vnNorSPosT[i]; j<=vnNorEPosT[i]; j++)
		{
			int nPos = j-nPivotSPos;
			if(nPos > 500)	continue;
			if(vsNorSeqT[i][j-vnNorSPosT[i]] == 'A')		pnNorA[nPos]++;
			if(vsNorSeqT[i][j-vnNorSPosT[i]] == 'C')		pnNorC[nPos]++;
			if(vsNorSeqT[i][j-vnNorSPosT[i]] == 'G')		pnNorG[nPos]++;
			if(vsNorSeqT[i][j-vnNorSPosT[i]] == 'T')		pnNorT[nPos]++;
			pnNorAll[nPos]++;
		}
	}

	for(int i=0; i<vnNorSPosT.size(); i++)
	{
		if(vnNorEPosT[i]-vnNorSPosT[i]+1 < nMidLen*0.8)		continue;
		int nCntMismatch = 0;
		for(int j=vnNorSPosT[i]; j<=vnNorEPosT[i]; j++)
		{
			if(j-vnNorSPosT[i] >= vsNorSeqT[i].size())	break;
			if(vsNorSeqT[i][j-vnNorSPosT[i]] == 'D')	continue;
			if(vsNorSeqT[i][j-vnNorSPosT[i]] != sRefSeq[j-nPivotSPos])		nCntMismatch++;
		}
		if(nCntMismatch > min(nMidMutation*2, 6))	continue;

		vnNorSPos.push_back(vnNorSPosT[i]);
		vnNorEPos.push_back(vnNorEPosT[i]);
		vsNorSeq.push_back(vsNorSeqT[i]);
	}
	for(int i=0; i<vnVarSPosT.size(); i++)
	{
	//	cout << vnVarSPosT[i] << "-" << vnVarEPosT[i] << "\t" << vsVarSeqT[i].size() << endl;	
		
		if(vnVarEPosT[i]-vnVarSPosT[i]+1 < nMidLen*0.8)		continue;
		int nCntMismatch = 0;
		for(int j=vnVarSPosT[i]; j<=vnVarEPosT[i]; j++)
		{
			if(j-vnVarSPosT[i] >= vsVarSeqT[i].size())	break;
			if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'D')	continue;
			
			int nPos = j-nPivotSPos;
			if(nPos > 500)		continue;
			if(pnVarAll[nPos] < 4)								continue;
			
			double dVarVaf = 0;
			double dNorVaf = 0;
			if(pnVarAll[nPos] > 0)
			{
				if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'A')		 dVarVaf = (double)pnVarA[nPos]/(double)pnVarAll[nPos];
				if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'C')		 dVarVaf = (double)pnVarC[nPos]/(double)pnVarAll[nPos];
				if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'G')		 dVarVaf = (double)pnVarG[nPos]/(double)pnVarAll[nPos];
				if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'T')		 dVarVaf = (double)pnVarT[nPos]/(double)pnVarAll[nPos];
			}
			if(pnNorAll[nPos] > 0)
			{
				if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'A')		 dNorVaf = (double)pnNorA[nPos]/(double)pnNorAll[nPos];
				if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'C')		 dNorVaf = (double)pnNorC[nPos]/(double)pnNorAll[nPos];
				if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'G')		 dNorVaf = (double)pnNorG[nPos]/(double)pnNorAll[nPos];
				if(vsVarSeqT[i][j-vnVarSPosT[i]] == 'T')		 dNorVaf = (double)pnNorT[nPos]/(double)pnNorAll[nPos];
			}

			if(dVarVaf >= 0.75 && dNorVaf <= 0.25)		continue;

			if(vsVarSeq[i][j-vnVarSPosT[i]] != sRefSeq[nPos])		nCntMismatch++;
		}
		//cout << nCntMismatch << endl;
		if(nCntMismatch > min(nMidMutation*2 + 1, 6))	continue;
		
		vnVarSPos.push_back(vnVarSPosT[i]);
		vnVarEPos.push_back(vnVarEPosT[i]);
		vsVarSeq.push_back(vsVarSeqT[i]);
	}

	// calculate read info.
	nSumNor = 0;
	nSumVar = 0;
	for(int i=0; i<vnNorSPosT.size(); i++)		nSumNor += vnNorEPosT[i]-vnNorSPosT[i]+1;
	for(int i=0; i<vnVarSPosT.size(); i++)		nSumVar += vnVarEPosT[i]-vnVarSPosT[i]+1;
	nCntNor = vnNorSPosT.size();
	nCntVar = vnVarSPosT.size();
	return true;
}
*/
/*
bool CRD::AlleleCount(bam_index_t *bamIndex, bamFile &finBam, 
		int nChr, int nPos, string sRef, string sAlt, int nRepeatCnt, int nRepeatLastLen, int &nPivotSPos, 
		vector<int> &vnNorSPos, vector<int> &vnNorEPos, vector<int> &vnVarSPos, vector<int> &vnVarEPos, 
		vector<string> &vsNorSeq, vector<string> &vsVarSeq)
{

	//CASE : seek start position of aligned reads in the bound
	bam_iter_t bamIter;
	bam1_t *b;
	b = bam_init1();
	bamIter = bam_iter_query(bamIndex, nChr, nPos-1, nPos);

	int nRet;

	int nCntVar = 0;
	int nCntNor = 0;
	while((nRet = bam_iter_read(finBam, bamIter, b)) >= 0)
	{
		//calc all depth
		int nSPos, nEPos;
		string sSeq;
		bool bIsVar = GetReadAlign(nSPos, nEPos, sSeq, nPos, sRef, sAlt, b);
		
		if(bIsVar)	nCntVar++;
		else		nCntNor++;
		

		//modify Pos depending Indel
		int nPosEdit = sRef.size()-sAlt.size();
		if(bIsVar)		nEPos -= nPosEdit;

		if(nSPos >= nPos-abs(nPosEdit))										continue;
//		if(nEPos <= nPos+abs(nPosEdit))										continue;
		if(nEPos <= nPos+abs(nPosEdit)*(nRepeatCnt)+nRepeatLastLen)			continue;


		if(nPivotSPos == -1 || (nPivotSPos != -1 && nPivotSPos > nSPos))	nPivotSPos = nSPos;		

		//calc variants depth
		if(bIsVar)
		{
			vnVarSPos.push_back(nSPos);
			vnVarEPos.push_back(nEPos);
			vsVarSeq.push_back(sSeq);
		}
		else
		{
			vnNorSPos.push_back(nSPos);
			vnNorEPos.push_back(nEPos);
			vsNorSeq.push_back(sSeq);
		}
	}

	bam_iter_destroy(bamIter);
	bam_destroy1(b);


	return true;
}


bool CRD::CalcPearsonCorr(int *pnAllDp, int *pnVarDp, double &dCorr, double &dPval)
{
	int nEPosT = 1000;
	for(int i=1000; i>=0; i--)
	{
		if(pnAllDp[i] != 0)	
		{
			nEPosT = i;
			break;
		}
//		cout << pnAllDp[i] << "\t";
	}
//	cout << endl;
//	for(int j=0; j<nEPosT; j++)		cout << pnVarDp[j] << "\t";
//	cout << endl;


	double *pdX = new double[nEPosT];
	double *pdY = new double[nEPosT];
	for(int i=0; i<nEPosT; i++)
	{
		pdX[i] = (double)pnAllDp[i];
		pdY[i] = (double)pnVarDp[i];
	}


	alglib::real_1d_array x;
	x.setcontent(nEPosT, pdX);
	alglib::real_1d_array y;
	y.setcontent(nEPosT, pdY);
	
	delete [] pdX;
	delete [] pdY;

	double dPearsonCorr = alglib::pearsoncorr2(x,y);
	alglib::ae_int_t nSize = nEPosT;
	double a, b, c;
	alglib::pearsoncorrelationsignificance(dPearsonCorr, nSize, a,b,c);
	

	dCorr = dPearsonCorr*dPearsonCorr;
	if(dPearsonCorr < 0)	dCorr = 0;
	if(dCorr < 0)	dCorr = 0;
	if(dCorr > 1)	dCorr = 1;
	dPval = a;
	return true;
}



bool CRD::GetReadAlign(int &nSPos, int &nEPos, string &sSeq, int nPos, string sRef, string sAlt, bam1_t *b)
{
	bool bIsVar = false;
	int nCigar[1000];
	char cSeq[1000];
	
	memcpy(nCigar, b->data+b->core.l_qname, b->core.n_cigar*4);
	memcpy(cSeq, b->data + b->core.l_qname + b->core.n_cigar*4, (b->core.l_qseq+1)/2);



	//resolve cigar
	int nRefPos = b->core.pos;
	int nSeqPos = 0;

	if(b->core.n_cigar == 1) 	// just one operation, save a loop
	{
		int op = nCigar[0] & 0x0000000f;
		int l = nCigar[0] >> 4;
		nSPos = nRefPos + 1;
		nEPos = nRefPos + l;
		
		if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
		{
			int nVarPosInSeq = nPos - (nRefPos+1);
			unsigned char base = bam1_seqi(cSeq, nVarPosInSeq); // base
			if(sRef.size() == 1 && sAlt.size() == 1)
			{			
				if((base == 1 && sAlt == "A") || (base == 2 && sAlt == "C"))	bIsVar = true;
				if((base == 4 && sAlt == "G") || (base == 8 && sAlt == "T"))	bIsVar = true;
			}
			
			//make tail seq
			for(int i=0; i<l; i++)		
			{	
				base = bam1_seqi(cSeq, i); // base

				if(base == 1)		sSeq.push_back('A');
				if(base == 2)		sSeq.push_back('C');
				if(base == 4)		sSeq.push_back('G');
				if(base == 8)		sSeq.push_back('T');
				if(base == 15)		sSeq.push_back('N');
			}
		}
	}
	else
	{
		nSPos = nRefPos+1;
		nEPos = nRefPos;	

		for(int i=0; i<b->core.n_cigar; i++)
		{ 
			int op = nCigar[i] & 0x0000000f;
			int l = nCigar[i] >> 4;
		
			if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
			{
				unsigned char base;
				int nVarPosInSeq = nSeqPos + (nPos - (nRefPos+1));
				
				if(nRefPos+1 <= nPos && nPos <= nRefPos+l)
				{
					base = bam1_seqi(cSeq, nVarPosInSeq); // base

					if(sRef.size() == 1 && sAlt.size() == 1)
					{
						if((base == 1 && sAlt == "A") || (base == 2 && sAlt == "C"))	bIsVar = true;
						if((base == 4 && sAlt == "G") || (base == 8 && sAlt == "T"))	bIsVar = true;
					}
				}
				
				//make tail seq
				for(int j=0; j<l; j++)		
				{	
					base = bam1_seqi(cSeq, nSeqPos+j);
					if(base == 1)		sSeq.push_back('A');
					if(base == 2)		sSeq.push_back('C');
					if(base == 4)		sSeq.push_back('G');
					if(base == 8)		sSeq.push_back('T');
					if(base == 15)		sSeq.push_back('N');
				}
				
				nRefPos += l;
				nSeqPos += l;
				nEPos += l;
			}
			else if(op == BAM_CDEL)
			{
				int nVarPosInSeq = nSeqPos + (nPos - (nRefPos+1));
				
				if(nRefPos+1 == nPos+1 && l == sRef.size()-sAlt.size())
				{
					unsigned char base = bam1_seqi(cSeq, nSeqPos-1); // base
					if((base == 1 && sAlt == "A") || (base == 2 && sAlt == "C"))	bIsVar = true;
					if((base == 4 && sAlt == "G") || (base == 8 && sAlt == "T"))	bIsVar = true;
				}
				
				for(int j=0; j<l; j++)		sSeq.push_back('D');
				
				nRefPos += l;
				nEPos += l;
			}
			else if(op == BAM_CREF_SKIP)
			{
			//	cout << "R" << l << endl;
				fprintf( stderr, "N" );
				fprintf( stderr, "(%d,%d): ", nRefPos, nSeqPos);
				nRefPos += l;
				nSeqPos += l;
				nEPos += l;
			}
			else if(op == BAM_CINS)
			{
			//	cout << "I" << l << endl;
				if(nRefPos+1 == nPos+1 && l == sAlt.size()-sRef.size())
				{
					//make insertion seq
					string sSeqTemp = "";
					for(int m=0; m<l; m++)
					{
						unsigned char base = bam1_seqi(cSeq, nSeqPos+m); //base
						switch(base)
						{
							case 1: sSeqTemp.push_back('A');    break;
							case 2: sSeqTemp.push_back('C');    break;
							case 4: sSeqTemp.push_back('G');    break;
							case 8: sSeqTemp.push_back('T');    break;
							case 15:sSeqTemp.push_back('N');    break;
							default:                break;
						}
					}

					if(sSeqTemp == sAlt.substr(1))
					{
						unsigned char base = bam1_seqi(cSeq, nSeqPos-1); // base
						if((base == 1 && sRef == "A") || (base == 2 && sRef == "C"))	bIsVar = true;
						if((base == 4 && sRef == "G") || (base == 8 && sRef == "T"))	bIsVar = true;
		
				//		cout << nPos << "\t" << nRefPos << "\t" << nSeqPos << "\t" << bIsVar << endl;
					}
				}
				nSeqPos += l;
			}
			else if(op == BAM_CSOFT_CLIP)
			{
				//cout << "S" << l << endl;
				nSeqPos += l;	
			}
			else if(op == BAM_CHARD_CLIP){}
		}
	}
	return bIsVar;
}




*/

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
















