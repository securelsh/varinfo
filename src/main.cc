 /**************************************************************************
  * *
  * *       Syntekabio, inc.
  * *
  * *  Copyright(c) 2018 Syntekabio, inc., 187, Techro 2-Gil, 
  * *  Yuseong-Gu, Daejeon, Korea.
  * *  All rights are reserved. No part of this work covered by the copyright
  * *  may be reproduced, stored in retrieval systems, in any form or by any
  * *  means, electronics, mechanical, photocopying, recording or otherwise,
  * *  without the prior permission of Syntekabio.
  * *
  * *  FILE NAME : main.cc
  * *  LAST MODIFIED DATE : 2018-12-29
  * *
  * *  AUTHOR : Sunho Lee
  * *  E-MAIL : shlee@syntekabio.com
  * *
  * *  DESCRIPTION : source file for VARinfo
  * *
  * **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <iostream>
#include <getopt.h>
#include <stdexcept>

//#include "../header/bam_to_fastq.h"
//#include "../header/gen_fasta.h"
//#include "../header/determine_seq.h"
#include "../header/varinfo.h"
//#include "../header/fisher.h"


using namespace std;

int parse_opt(const int argc, char **argv, int &nStatus, string &sBamFile, 
	string &sInputFile, string &sOutputFile, bool &bIsMod, int &nCntThread, bool &bIsDebug);
bool manual()
{
	cout << "====================================== RDscan manual ====================================" << endl;
	cout << "1. COMMANDs"											<< endl;
	cout << "  (1) Output Variant Info" 									<< endl;
    cout << "      ./varinfo -b ${Bamfile} -i ${Adifile} -o ${Output}" 	<< endl;

	cout << endl;
	cout << "2. OPTIONs" << endl;
	cout << "  (1) # of Thread                     -t    [unsigned int], default=1"				<< endl;
	cout << "  (2) consecurive SNPs to Indel       -m"	<< endl;
	cout << "==========================================================================================" << endl << endl;

	return true;
}


int main(int argc, char* argv[])
{
	int nStatus = 1;			// 1:vcf 
	string sBamFile = "";			//bam file
	string sInputFile;
	string sOutputFile;
	int nCntThread = 1;
	bool bIsDebug = false;
	bool bIsMod = false;


	try
	{
		parse_opt(argc, argv, nStatus, sBamFile, sInputFile, sOutputFile, bIsMod, nCntThread, bIsDebug); 
		if(nStatus == 0)			manual();
		else
		{
			//construct class & init
			CINFO INFO(sBamFile, sInputFile, sOutputFile, bIsMod, nCntThread, bIsDebug);

		
			INFO.PrintCommonInfo();
			cout << endl;
 
			//read files (ref, input)
			INFO.ReadInput();
		
			if(bIsMod)		INFO.ModVariant();

			//calculate read distribution
			INFO.CalcInfo();

			//report
			//RD.Report();
		}
	}
	catch (std::logic_error const& s)
	{
		cout << s.what();
		exit(EXIT_FAILURE);
	}
	catch (bad_alloc &ba)
	{
		cout << ba.what() << endl;
		exit(EXIT_FAILURE);
	}
	catch (...)
	{
		cout << "ERROR: Unexpected Exception" << endl;
		exit(EXIT_FAILURE);
	}


	return false;
}


int parse_opt(const int argc, char **argv, int &nStatus, string &sBamFile, 
string &sInputFile, string &sOutputFile, bool &bIsMod, int &nCntThread, bool &bIsDebug)
{
	bool bBOption = false;
	bool bIOption = false;
	bool bOOption = false;

	struct option long_options[] =
	{
		{"bamFile", 		1, 0, 'b'},		//0
		{"input",			1, 0, 'i'},		//3
		{"output",			1, 0, 'o'},
		{"modigy",			1, 0, 'm'},		//4	
		{"num_thread", 		1, 0, 't'},		//8
		{0, 0, 0, 0}
	};

	int c;
	
	while (1) 
	{   
		int option_index = 0;
		c = getopt_long(argc, argv, "b:i:o:t:mq", long_options, &option_index);
		if (c == -1)    break;
		char *pcEnd;
		switch (c) 
		{   
			case 'b':
				bBOption = true;
				sBamFile = optarg;
				break;
			case 'm':
				bIsMod = true;
				break;
			case 'i':
				bIOption = true;
				sInputFile = optarg;
				break;
			case 'o':
				bOOption = true;
				sOutputFile = optarg;
				break;
			case 't':
				nCntThread = strtol(optarg, &pcEnd, 10);
				break;
			case 'q':
				bIsDebug = true;
				break;
		}
	}

	if(!bBOption || !bIOption || !bOOption)		nStatus = 0;
	if(nCntThread<1) nCntThread = 1;
	return true;
}










