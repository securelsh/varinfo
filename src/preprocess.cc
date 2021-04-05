#include "../header/varinfo.h"
#include "../header/bio_files.h"
#include "../header/statisticaltest.h"


using namespace std;
pthread_mutex_t mutex_temp = PTHREAD_MUTEX_INITIALIZER;

/*

- To Do:
  -- get read info for multiple loci
  
- Changelog:
  -- Date Mar-30-2021 SeokCholHong shulkhorn@gmail.com
     --- GetRead

*/

bool CINFO::ModVariant()
{
	for(int i=0; i<this->m_Input.vbIsVcf.size()-1; i++)
	{
		int nIdx = i;
		while(this->m_Input.vnPos[nIdx] + 1 == this->m_Input.vnPos[nIdx+1]
		&& this->m_Input.vsRef[nIdx].size() == 1 && this->m_Input.vsAlt[nIdx].size() == 1
		&& this->m_Input.vsRef[nIdx+1].size() == 1 && this->m_Input.vsAlt[nIdx+1].size() == 1)
		{
			nIdx++;
			if(nIdx == this->m_Input.vbIsVcf.size()-1)	break;
		}


		if(i != nIdx)			
		{
			this->DelIns(i, nIdx);
			i = nIdx;

		}


	}
	return false;
}


bool CINFO::DelIns(int nSIdx, int nEIdx)
{
	bamFile finBam;
	finBam = bam_open(m_sBamFile.c_str(), "r");
	if(finBam == NULL)	throw	std::logic_error("ERROR: Can not open .bam file");

	bam_header_t *bamHeader;
	bamHeader = bam_header_read(finBam);
	if(!bamHeader)		throw std::logic_error("ERROR: Fail to read the header of a bam file");

	bam_index_t *bamIndex;
	bamIndex = bam_index_load(m_sBamFile.c_str());
	if(!bamIndex)		throw std::logic_error("ERROR: Fail to read .bai file");

	string sChr = this->m_Input.vsChr[nSIdx];
	int nChr = ConvertChrToTid(sChr, bamHeader);
	int nPos = this->m_Input.vnPos[nSIdx];

	if(nChr == -1)
	{
		cout << sChr << "(" << nChr << "):" << nPos << endl;
		throw std::logic_error("Not available Chr");
	}



	for(int i=nSIdx; i<nEIdx; i++)
	{
		nPos = this->m_Input.vnPos[i];
		string sVarSeq = this->m_Input.vsAlt[i];
		sVarSeq += this->m_Input.vsAlt[i+1];
		// extract reads
		bam_iter_t bamIter;
		bam1_t *b;
		b = bam_init1();
		bamIter = bam_iter_query(bamIndex, nChr, nPos-1, nPos);
		int nRet;
		int nOO = 0;
		int nAO = 0;
		int nOA = 0;
		int nAA = 0;
		while((nRet = bam_iter_read(finBam, bamIter, b)) >= 0)
		{
			string sSeq = "";
			for(int j=nPos; j<=nPos+1; j++){
				int nLoc = j - (b->core.pos+1);
				if(0<=nLoc){
					if(GetRead(b, nLoc)=="") break;
					sSeq+=GetRead(b, nLoc);
				}
			}
			//cout << sSeq << " ";
			if(sSeq[0] == sVarSeq[0] && sSeq[1] == sVarSeq[1])		nAA++;
			if(sSeq[0] == sVarSeq[0] && sSeq[1] != sVarSeq[1])		nAO++;
			if(sSeq[0] != sVarSeq[0] && sSeq[1] == sVarSeq[1])		nOA++;
			if(sSeq[0] != sVarSeq[0] && sSeq[1] != sVarSeq[1])		nOO++;


		}
		int nOriOO = nOO;
		int nOriOA = nOA;
		int nOriAO = nAO;
		int nOriAA = nAA;
		//cout << endl;
		if(m_bIsDebug)	cout << "(" << nOO << "," << nAO << "," << nOA << "," << nAA << ")" << endl;
		if(m_bIsDebug)	cout << m_Input.vdVaf[i] << "\t" << m_Input.vdVaf[i+1] << endl;
		nAO *= m_Input.vdVaf[i];
		nOA *= m_Input.vdVaf[i+1];
		nOO *= m_Input.vdVaf[i] * m_Input.vdVaf[i+1];
		double dPvalue = (double)STATTEST::GetFisherPvalue(nOO, nAO, nOA, nAA);
		if(m_bIsDebug)	cout << "P-value: " << dPvalue << " (" << nOO << "," << nAO << "," << nOA << "," << nAA << ")" << endl;

		// determine whether SNPs or Del-Ins
		if(dPvalue < 0.05)
		{
			string sRaw = "chr " + this->m_Input.vsChr[i] + "_" + this->to_string(this->m_Input.vnPos[i]);
			sRaw += "\t" + this->to_string(nOriOO) + " " + this->to_string(nOriAA);
			sRaw += "\t" + this->m_Input.vsRef[i] + this->m_Input.vsRef[i+1] + " " + sVarSeq;

			this->m_Input.vbIsPass[i] = false;
			this->m_Input.vnPos[i+1]--;
			this->m_Input.vsRef[i+1] = this->m_Input.vsRef[i] + this->m_Input.vsRef[i+1];
			this->m_Input.vsAlt[i+1] = this->m_Input.vsAlt[i] + this->m_Input.vsAlt[i+1];
			this->m_Input.vdVaf[i+1] = (double)nOriAA / (double)(nOriOO+nOriAA);
			

			char *pcEnd;
			vector<string> vsWord1, vsWord2, vsT1, vsT2;
			Parsing(vsWord1, m_Input.vsRaw[i], "\t");
			Parsing(vsWord2, m_Input.vsRaw[i+1], "\t");
			Parsing(vsT1, vsWord1[2], " ");
			Parsing(vsT2, vsWord2[2], " ");

			if(vsT1[3] == "2" || vsT2[3] == "2")	sRaw += " | 2 | - | ";
			else									sRaw += " | 3 | - | ";
			double dScore = (strtod(vsT1[7].c_str(), &pcEnd) + strtod(vsT2[7].c_str(), &pcEnd)) / 2.0;
			sRaw += this->to_string(dScore) + " | - 0 | -";

			cout << sRaw << endl;
			this->m_Input.vsRaw[i+1] = sRaw;	
		}
	}


	bam_close(finBam);

	return false;
}

string CINFO::GetRead(bam1_t *b, int nPos)
{
	string sSeq="";
	//cigar
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
				if(nPos<nPreAlgn)	break;
			}
			else if(op == BAM_CINS){
				if(nPos<nPreAlgn)	break;
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CDEL){
				if(nPreAlgn<=nPos && nPos<nPreAlgn+nAlgn){
					return sSeq;
				}
				else if(nPos>=nPreAlgn+nAlgn){
					nPos -= nAlgn;
				}
				else if(nPos<nPreAlgn)	break;
			}
			else if(op == BAM_CSOFT_CLIP || op == BAM_CPAD){
				if(nPos<nPreAlgn)	break;
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CREF_SKIP){
				if(nPos<nPreAlgn)	break;
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
	if(nPos>=b->core.l_qseq)	return sSeq;

	unsigned char ucBase = bam1_seqi(bam1_seq(b), nPos);
	if     (ucBase == 0x01)	sSeq="A";
	else if(ucBase == 0x02)	sSeq="C";
	else if(ucBase == 0x04)	sSeq="G";
	else if(ucBase == 0x08)	sSeq="T";
	else if(ucBase == 0x0F)	sSeq="N";
	else
		throw std::logic_error("ERROR: sequence error. It deosn't belong to ACGTN (in function GetRead())");

	return sSeq;
}



CINFO::CINFO(string sBamFile, string sInputFile, string sOutputFile, bool bIsMod, int nCntThread, double dVaf, bool bIsDebug)
{
	m_sBamFile = sBamFile;
	m_sInputFile = sInputFile;
	m_sOutputFile = sOutputFile;
	m_bIsMod = bIsMod;
	m_nCntThread = nCntThread;
	m_dVaf = dVaf;
	m_bIsDebug = bIsDebug;
}


bool CINFO::ReadInput()
{
	//read vcf file
	this->SetStartTime();
	this->ViewStatus(0,1,"Read Input File");
	
	ifstream fin;
	fin.open(m_sInputFile.c_str());
	if(!fin.is_open())		throw std::logic_error("Input file does not exists");

	string sLine;
	char *pcEnd;
	vector<string> vsWord;
	while(getline(fin, sLine) && sLine.size() > 0)
	{
		if(sLine[0] == '#')			this->m_Input.vsHeader.push_back(sLine);
		else
		{
			Parsing(vsWord, sLine, "\t");
			vector<string> vsDepth;
			Parsing(vsDepth, vsWord[1], " ");

			int nRefDp = strtol(vsDepth[0].c_str(), &pcEnd, 10);
			int nAltDp = strtol(vsDepth[1].c_str(), &pcEnd, 10);

			double dVaf = (double)nAltDp/(double)(nRefDp+nAltDp);
			if(dVaf < this->m_dVaf)	continue;

			this->m_Input.vdVaf.push_back(dVaf);

			if(vsWord[0].substr(0,4) == "chr ")			//Adiscan format
			{
				this->m_Input.vbIsVcf.push_back(false);
				vector<string> vsSplit;
				Parsing(vsSplit, vsWord[0].substr(4), "_");
				this->m_Input.vsChr.push_back(vsSplit[0]);
				this->m_Input.vnPos.push_back(strtol(vsSplit[1].c_str(), &pcEnd, 10));
				
				Parsing(vsSplit, vsWord[2], " ");
				this->m_Input.vsRef.push_back(vsSplit[0]);
				this->m_Input.vsAlt.push_back(vsSplit[1]);
			}
			else										//Vcf format
			{
				this->m_Input.vbIsVcf.push_back(true);
				this->m_Input.vsChr.push_back(vsWord[0]);
				this->m_Input.vnPos.push_back(strtol(vsWord[1].c_str(), &pcEnd, 10));
				this->m_Input.vsRef.push_back(vsWord[3]);
				this->m_Input.vsAlt.push_back(vsWord[4]);
	
			}
			this->m_Input.vsRaw.push_back(sLine);
			this->m_Input.vbIsPass.push_back(true);
		}
	}
	
	this->ViewStatus(1,1,"Read Input File");
	cout << endl;


	return false;
}




