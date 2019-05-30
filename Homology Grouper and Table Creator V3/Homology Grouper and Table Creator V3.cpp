// Homology Grouper and Table Creator V3.cpp : Defines the entry point for the console application.
// Author: Mohamed El-Walid - mze99d@mail.missouri.edu
// University of Missouri - Columbia
// Kathleen Newton Lab

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>     
#include <time.h> 

using namespace std;


int initUserInput(string*, string*, string*, string*, string*, int*, int*, int*);
int nucleoCount(string);
int countGenomes(string*, int*);

int processInputHomologies(string, int*, int***, int, int, string);
int countHomologies(string*, int*);
int readInData(string, int, int***);
int** removeRedundancies(int***, int*, int*, int);

//void groupingHomologies(int***, int, int, int);
//void generateRepresentation(int**, int***, int, int, int);
//void generateGroupings(int**, int, int);
//int seqPresent(int**, int, int, int);
//int checkNext(int**, int, int, int);
//bool doMatch(int**, int, int, int, int);

void quickSortArray(int**, int, int);
int partition(int**, int, int);
void arraySwap(int**, int, int);
int medianOfThreeIndex(int**, int, int, int);
void bubbleSortArray(int**, int, int);

void generateGroups(int**, int, int, string, int);
bool isInGroup(int**, int, int);
void finalGroups(string*, int);

int regionsExtraction(string, string);


const int ALLOWEDGAP = 5000;
//try varius mitogaps; smaller prob better; 5bp maybe
const int ALLOWEDMITOGAP = 100;

const int CONDITIONALGAP = 20000;
const int CONDITIONALMITOGAP = 100;



int main()
{
	int sourceSize, numOfGenomes, numOfHomologies, numOfRegions = 0;
	srand(time(NULL));

	//establish file names ------------------------------------------------------------------------------------------------------------------------------------------------
	string  regionOutput, dataOutput, listFile, sourceFileName, sourceName;

	if(initUserInput(&regionOutput, &dataOutput, &listFile, &sourceFileName, &sourceName, &sourceSize, &numOfGenomes, &numOfHomologies) == 1){
		cout << "initUserInput function couldn't open something.";
		cout << "\nEnter any key to continue....";
		getchar();
		getchar();
		return 0;
	}

	cout << "Scaffold Size is: " << sourceSize << endl;
	cout << listFile << endl << numOfHomologies << endl << numOfGenomes;
	getchar();

	//proccess input data --------------------------------------------------------------------------------------------------------------------------------------------------------
	
	int **dataArray = NULL;
	if (processInputHomologies(listFile, &numOfHomologies, &dataArray, numOfGenomes, sourceSize, sourceFileName) == 1)
	{
		cout << "processInputHomologies died";
		return 1;
	}

	cout << "out of processInput";

	for (int i = 0; i < 6; i++)
		delete[] dataArray[i];
	delete dataArray;
	
    return 0;
}

int initUserInput(string *regionOutput, string *dataOutput, string *listFile, string *scaffoldFileName, string *scaffoldName, int *scaffoldSize, int *numOfGenomes, int *numOfHomologies)
{
	cout << "This Program requires a source sequence file name and text list of individual genome containing files.\n\nEnter the scaffold name: ";

	cin >> *scaffoldName;

	cout << "Input the name of your list file (inlcude file type): ";
	cin >> *listFile;

	*regionOutput = *scaffoldName + "_Scaffold_Homologous_Groups.txt";
	*dataOutput = *scaffoldName + "_Genome_Homo_Grouping.txt";

	cout << "Scaffold File Name (inlcude file type): ";
	cin >> *scaffoldFileName;

	*scaffoldSize = nucleoCount(*scaffoldFileName);
	
	if (*scaffoldSize == -1)
		return 1;

	if (countGenomes(listFile, numOfGenomes ) == -1)
		return 1;

	if (countHomologies(listFile, numOfHomologies) == -1)
		return 1;

	return 0;
}

int nucleoCount(string fileName)
{
	ifstream scaffold;
	bool inSeq = false;
	int count = 0;
	char hold;

	scaffold.open(fileName);

	if (!scaffold.is_open())
	{
		cout << "Scaffold file could not be opened. Returning error 1.\n";
		return -1;
	}

	while (!scaffold.eof())
	{
		scaffold.get(hold);
		tolower(hold);

		if (((hold > 64 && hold < 90) || (hold > 96 && hold < 123)) && inSeq == true)
			count++;

		if (hold == '\n' || hold == 'j')
			inSeq = true;
	}

	if (scaffold.is_open())
		scaffold.close();
	return count;
}

int countGenomes(string *listFile, int *numOfGenomes) {
	ifstream input(*listFile);
	*numOfGenomes = 0;

	if (!input.is_open())
		return -1;
	cout << "numOfGenomes: " << *numOfGenomes << endl;;
	string hold;
	while (!input.eof()){
		input >> hold;
		*numOfGenomes = *numOfGenomes + 1;
		cout << hold << " " << *numOfGenomes <<" ";
	}
	
	input.close();

	return 0;
}

int countHomologies(string *listFile, int *homologiesCount)
{
	int mStart, mStop, sStart, sStop;
	string currGenome, currFileName;
	*homologiesCount = 0;

	ifstream listInput(*listFile);
	ifstream dataInput;
	if (!listInput.is_open())
		return 1;

	while (!listInput.eof())
	{
		listInput >> currGenome;
		currFileName = currGenome + ".txt";
		dataInput.open(currFileName);

		if (dataInput.is_open()){
			while (!dataInput.eof()){
				dataInput >> mStart >> mStop >> sStart >> sStop;
				if(mStart >= 0 && mStop >= 0 && sStart >= 0 && sStop >= 0)
					*homologiesCount = *homologiesCount + 1;
			}
		}
		else
			cout << "File \"" + currFileName + "\" could not be opened" << endl;

		dataInput.close();
	}

	if(listInput.is_open())
		listInput.close();

	return 0;
}

int processInputHomologies(string listFile, int *numOfHomologies, int ***array, int numOfGenomes, int sizeOfSource, string sourceFileName)
{
	*array = new int*[6];
	int **dataArray = *array;

	for (int i = 0; i < 6; i++)
		dataArray[i] = new int[*numOfHomologies];

	int *homologiesPerGenome = new int[numOfGenomes];

	if (readInData(listFile, *numOfHomologies, &dataArray) == 1)
	{
		for (int i = 0; i < 6; i++)
			delete[] dataArray[i];
		delete dataArray;
		return 1;
	}

	cout << dataArray << " " << dataArray[3][0] << " " << *numOfHomologies << endl;
	int **newArray = removeRedundancies(&dataArray, homologiesPerGenome, numOfHomologies, numOfGenomes);
	cout << newArray << " " << newArray[3][0] << " " << *numOfHomologies << endl;
	
	ofstream output("output.txt");
	for (int i = 0; i < *numOfHomologies; i++) {
		output << newArray[0][i] << " | " << newArray[1][i] << " | " << newArray[2][i] << " | " << newArray[3][i] << " | " << newArray[4][i] << endl;
	}
	output.close();

	int currGenomeStart = 0;
	for (int i = 0; i < numOfGenomes; i++)
	{
		quickSortArray(newArray, currGenomeStart, (currGenomeStart + homologiesPerGenome[i]) - 1);
		bubbleSortArray(newArray, currGenomeStart, (currGenomeStart + homologiesPerGenome[i]) - 1);
		currGenomeStart += homologiesPerGenome[i];
	}

	output.open("sortedOutput.txt");
	for (int i = 0; i < *numOfHomologies; i++) {
		output << newArray[0][i] << " | " << newArray[1][i] << " | " << newArray[2][i] << " | " << newArray[3][i] << " | " << newArray[4][i] << endl;
	}
	output.close();

	string fileName = "groupings.txt";
	output.open(fileName, ofstream::out | ofstream::trunc);
	output.close();

	currGenomeStart = 0;
	for (int i = 0; i < numOfGenomes; i++)
	{
		generateGroups(newArray, currGenomeStart, ((currGenomeStart + homologiesPerGenome[i]) - 1), fileName, sizeOfSource);
		cout << endl << currGenomeStart << " - " << (currGenomeStart + homologiesPerGenome[i] - 1);
		currGenomeStart += homologiesPerGenome[i];
	}

	
	finalGroups(&fileName, sizeOfSource);


	for (int i = 0; i < 5; i++)
		delete[] newArray[i];
	delete newArray;

	regionsExtraction(sourceFileName, fileName);

	return 0;
}

int readInData(string listFile, int numOfHomology, int ***array)
{
	int **dataArray = *array;
	int temp, counter = 0;
	int genomeCount = 0;
	int mStart, mStop, sStart, sStop;
	ifstream genomesFile;
	ifstream dataFile;
	string genomeHold;

	genomesFile.open(listFile);

	if (!genomesFile.is_open())
		return 1;

	while (!genomesFile.eof())
	{
		genomesFile >> genomeHold;
		dataFile.open(genomeHold + ".txt");
		if (dataFile.is_open())
		{
			while (!dataFile.eof())
			{
				dataFile >> mStart >> mStop >> sStart >> sStop;
				if (mStart > 0 && mStop > 0 && sStart > 0 && sStop > 0) {
					dataArray[0][counter] = genomeCount;
					dataArray[1][counter] = mStart;
					dataArray[2][counter] = mStop;
					dataArray[3][counter] = sStart;
					dataArray[4][counter] = sStop;
					dataArray[5][counter] = 0;

					if (dataArray[3][counter] > dataArray[4][counter])
					{
						temp = dataArray[3][counter];
						dataArray[3][counter] = dataArray[4][counter];
						dataArray[4][counter] = temp;

						temp = dataArray[1][counter];
						dataArray[1][counter] = dataArray[2][counter];
						dataArray[2][counter] = temp;
					}
					counter++;
				}
			}
		}
		else
			return 1;

		dataFile.close();
		genomeCount++;
	}
	genomesFile.close();

	ofstream test;
	test.open("tester.txt");

	for (int i = 0; i < numOfHomology; i++)
		test << dataArray[0][i] << " " << dataArray[1][i] << " " << dataArray[2][i] << " " << dataArray[3][i] << " " << dataArray[4][i] << " " << dataArray[5][i] << "\n";

	test.close();

	return 0;
}

int** removeRedundancies(int ***array, int *homologiesPerGenome, int *numOfHomologies, int numOfGenomes)
{
	int genomeNum = 0;
	int colStop, colStart = 0;
	int indexHold, sourceStart, sourceStop, sourceStartCompare, sourceStopCompare;
	int **dataArray = *array;
	int redundantCounter = 0;
	

	for (int i = 0; i < numOfGenomes; i++)
		homologiesPerGenome[i] = 0;
	cout << "into while loop\n";
	while (genomeNum < numOfGenomes)
	{
		for (int i = colStart; dataArray[0][i] == genomeNum && i < *numOfHomologies; i++)
		{
			colStop = i + 1;
			sourceStart = dataArray[3][i];
			sourceStop = dataArray[4][i];
			indexHold = i;

			for (int j = i + 1; dataArray[0][j] == genomeNum && j < *numOfHomologies; j++)
			{
				if (dataArray[5][j] == 0 && dataArray[0][j] == dataArray[0][indexHold])
				{
					sourceStartCompare = dataArray[3][j];
					sourceStopCompare = dataArray[4][j];

					if ((sourceStart <= sourceStartCompare && sourceStop >= sourceStopCompare) ||
						(sourceStartCompare <= sourceStart && sourceStopCompare >= sourceStop)){
						if (sourceStop - sourceStart < sourceStopCompare - sourceStartCompare){
							dataArray[5][indexHold] = 1;
							indexHold = j;
							sourceStop = sourceStopCompare;
							sourceStart = sourceStartCompare;
							redundantCounter++;
						}
						else
						{
							dataArray[5][j] = 1;
							redundantCounter++;
						}
					}
				}
			}
		}

		genomeNum++;
		colStart = colStop;
	}
	cout << "out of while loop\n";

	int newNumOfHomologies = *numOfHomologies - redundantCounter;

	int **newDataArray = new int*[5];
	int j = 0;

	for (int i = 0; i < 5; i++)
		newDataArray[i] = new int[newNumOfHomologies];
	
	for (int i = 0; i < *numOfHomologies; i++)
	{
		if (dataArray[5][i] == 0)
		{
			cout << " " << j << " ";
			if (j < newNumOfHomologies) {
				newDataArray[0][j] = dataArray[0][i];
				newDataArray[1][j] = dataArray[1][i];
				newDataArray[2][j] = dataArray[2][i];
				newDataArray[3][j] = dataArray[3][i];
				newDataArray[4][j] = dataArray[4][i];
				(homologiesPerGenome[newDataArray[0][j]]) += 1;
				j++;
			}
		}
	}

	for (int i = 0; i < 6; i++) {
		cout << " " << dataArray[i] << " ";
	}

	cout << dataArray << endl;

	ofstream test;
	test.open("tester2.txt");
	for (int i = 0; i < newNumOfHomologies; i++) {
		test << newDataArray[0][i] << " | " << newDataArray[1][i] << " | " << newDataArray[2][i] << " | " << newDataArray[3][i] << " | " << newDataArray[4][i] << endl;
	}

	test.close();

	cout << "issue could be in delete";
	for (int i = 0; i < 6; i++) {
		if(dataArray[i])
			delete[] dataArray[i];
	}

	if(dataArray)
		delete dataArray;

	cout << "got past delete";
	*numOfHomologies = newNumOfHomologies;
	cout << "At end of Remove Redundancies";
	return newDataArray;
}

/*void groupingHomologies(int ***array, int numOfHomologies, int numOfGenomes, int sizeOfSource)
{
	int **dataArray = *array;

	//int **sourceRepresentation = new int*[numOfGenomes];

	//for (int i = 0; i < numOfGenomes; i++)
		//sourceRepresentation[i] = new int[sizeOfSource]();

	cout << "Going into generateRepresentation\n";
	//generateRepresentation(dataArray, &sourceRepresentation, numOfHomologies, numOfGenomes, sizeOfSource);

	//cout << "\n\nOut of representation\n";
	//cout << "Going into generateGroupings\n";
	//generateGroupings(sourceRepresentation, numOfGenomes, sizeOfSource);

	//for (int i = 0; i < numOfGenomes; i++)
		//delete[] sourceRepresentation[i];
	//delete sourceRepresentation;

	return;
}*/

/*void generateRepresentation(int **dataArray, int ***array, int numOfHomologies, int numOfGenomes, int sizeOfSource)
{
	
	int **sourceRepresentation = new int*[numOfGenomes];

	for (int i = 0; i < numOfGenomes; i++)
		sourceRepresentation[i] = new int[50];
	cout << "done Allocating" << endl;

	

	cout << setw(15) << dataArray[3][0] << setw(15) << *array << setw(15) << sourceRepresentation << setw(15) << numOfHomologies << setw(15) << numOfGenomes << setw(15) << sizeOfSource << endl;

	for (int i = 0; i < numOfHomologies; i++)
	{
		if ((dataArray[3][i] > 0) && (dataArray[4][i] < sizeOfSource)) 
		{
			for (int j = (dataArray[3][i]) + 1; j <= (dataArray[4][i]) - 3; j++) {
				sourceRepresentation[(dataArray[0][i]) - 1][j] = -1;
			}

			if ((sourceRepresentation[dataArray[0][i]][(dataArray[3][i]) - 1] == 0) && (sourceRepresentation[dataArray[0][i]][(dataArray[3][i])] == 0)) {
				sourceRepresentation[dataArray[0][i]][(dataArray[3][i]) - 1] = dataArray[1][i];
				sourceRepresentation[dataArray[0][i]][(dataArray[3][i])] = dataArray[2][i];
			}

			if ((sourceRepresentation[dataArray[0][i]][(dataArray[4][i]) - 1] == 0) && (sourceRepresentation[dataArray[0][i]][(dataArray[4][i]) - 2] == 0)) {
				sourceRepresentation[dataArray[0][i]][(dataArray[4][i]) - 1] = dataArray[2][i];
				sourceRepresentation[dataArray[0][i]][(dataArray[4][i]) - 2] = dataArray[1][i];
			}
		}
		cout << i << " ";
	}
	return;
}*/

/*void generateGroupings(int **sourceRepresentation, int numOfGenomes, int sizeOfSource)
{
	int startOne, stopOne, startTwo, stopTwo;
	int count;
	bool inRegion = false;
	int holdStart, holdStop;
	ofstream output;

	output.open("Groupings.txt");

	for (int i = 0; i < sizeOfSource; i++)
	{
		if (inRegion == false && seqPresent(sourceRepresentation, numOfGenomes, sizeOfSource, i) > 0)
		{
			holdStart = i;
			inRegion = true;
		}

		if (inRegion == true && seqPresent(sourceRepresentation, numOfGenomes, sizeOfSource, i) > 0 &&
			seqPresent(sourceRepresentation, numOfGenomes, sizeOfSource, i + 1) == 0)
		{
			holdStop = i;
			int hold = checkNext(sourceRepresentation, numOfGenomes, sizeOfSource, i + 1);
			
			if (hold > i + 1)
				i = hold;
			else
			{
				inRegion = false;
				output << setw(15) << holdStart << holdStop << endl;
				i = hold;
			}
		}

	}
	output.close();
	return;
}*/

/*int seqPresent(int **sourceRepresentation, int numOfGenomes, int sizeOfSource, int index)
{
	if (index >= sizeOfSource)
		return -1;

	int result = 0;

	for (int i = 0; i < numOfGenomes; i++)
	{
		if (sourceRepresentation[i][index] > 0)
			result = 1;
	}

	for (int i = 0; i < numOfGenomes; i++)
	{
		if (sourceRepresentation[i][index] == -1)
			result = 2;
	}

	return result;
}*/

/*bool doMatch(int **sourceRepresentation, int numOfGenomes, int sizeOfSource, int initial, int compared)
{
	int startOne, stopOne, startTwo, stopTwo;
	int counter = 0;
	bool result = false;

	if (compared >= sizeOfSource)
		return false;

	while (counter < numOfGenomes)
	{
		startOne = sourceRepresentation[counter][initial - 1];
		stopOne = sourceRepresentation[counter][initial];
		startTwo = sourceRepresentation[counter][compared];
		stopTwo = sourceRepresentation[counter][compared + 1];

		if (startOne > 0 && startTwo > 0 && stopOne > 0 && stopTwo > 0)
		{
			if ((stopOne > startTwo) && (stopOne < stopTwo))
			{
				return true;
			}
			else if ((startTwo > stopOne) && (startTwo < startOne))
				return true;
			else
			{
				int gap = startTwo - stopOne;
				if (gap < 0)
					gap = stopOne - startTwo;

				if (gap < 500 && gap >= 0)
					return true;
			}
		}
		counter++;
	}

	if (compared - initial < ALLOWEDGAP)
		return true;

	return false;
}*/

void quickSortArray(int **dataArray, int start, int stop) {
	//sort using dataArray[3][]
	
	if (start < stop)
	{	
		int div = partition(dataArray, start, stop);
		
		quickSortArray(dataArray, start, div - 1);  // Before pi
		quickSortArray(dataArray, div + 1, stop); // After pi
	}
}

int partition(int **dataArray, int start, int stop)
{
	if (stop == start + 1)
	{
		if (dataArray[3][start] > dataArray[3][stop])
			arraySwap(dataArray, start, stop);
		return start;
	}


	int left = start;
	int right = stop-1;

	int one, two, three, pivotValue;
	one = (rand() % (stop - left) + left) ;
	two = (rand() % (stop - left) + left);
	three = (rand() % (stop - left) + left);
	int pivotIndex = medianOfThreeIndex(dataArray, one, two, three);
	pivotValue = dataArray[3][pivotIndex];
	
	arraySwap(dataArray, pivotIndex, stop);

	while (left < right)
	{
		while (dataArray[3][left] < pivotValue)
			left++;

		while (dataArray[3][right] > pivotValue)
			right--;

		if (left < right)
		{
			arraySwap(dataArray, left, right);
			left++;
			right--;
		}
	}

	arraySwap(dataArray, left, stop);

	return left;
}

int medianOfThreeIndex(int **dataArray, int a, int b, int c)
{
	int temp;

	if (dataArray[3][a] > dataArray[3][b])
	{
		if (dataArray[3][a] > dataArray[3][c])
		{
			if (dataArray[3][b] > dataArray[3][c])
				return b;
			else
				return c;
		} 
		else 
			return a;
	} 
	else {
		if (dataArray[3][a] > dataArray[3][c])
			return a;
		else
		{
			if (dataArray[3][b] > dataArray[3][c])
				return c;
			else
				return b;
		}
	}

	return b;
}

void arraySwap(int **dataArray, int a, int b)
{
	if (a == b)
		return;

	int zero, one, two, three, four;

	zero = dataArray[0][a];
	one = dataArray[1][a];
	two = dataArray[2][a];
	three = dataArray[3][a];
	four = dataArray[4][a];
	

	dataArray[0][a] = dataArray[0][b];
	dataArray[1][a] = dataArray[1][b];
	dataArray[2][a] = dataArray[2][b];
	dataArray[3][a] = dataArray[3][b];
	dataArray[4][a] = dataArray[4][b];
	

	dataArray[0][b] = zero;
	dataArray[1][b] = one;
	dataArray[2][b] = two;
	dataArray[3][b] = three;
	dataArray[4][b] = four;
	

	return;
}

void bubbleSortArray(int **dataArray, int start, int stop) 
{
	bool isSorted;

	do {
		isSorted = true;
		for (int i = start; i < stop - 1; i++)
		{
			if (dataArray[3][i] > dataArray[3][i + 1])
			{
				arraySwap(dataArray, i, i + 1);
				isSorted = false;
			}
		}
	} while (isSorted == false);
}

void generateGroups(int **dataArray, int start, int stop, string fileName, int sizeOfSource)
{
	bool inRegion = false;
	int holdStart, holdStop;
	ofstream output;
	output.open(fileName, ios::out | ios::app);

	for (int i = start; i <= stop; i++)
	{
		if (inRegion == false)
		{
			holdStart = i;
			holdStop = i;
			inRegion = true;
		}
		else if (i == stop) {
			output << dataArray[3][holdStart] << setw(15) << dataArray[4][holdStop] << endl;
		}
		else
		{
			if (isInGroup(dataArray, holdStop, i))
				holdStop = i;
			else
			{
				output << dataArray[3][holdStart] << setw(15) << dataArray[4][holdStop] << endl;
				holdStart = holdStop = i;
			}
		}
		
	}

	output.close();
	return;
}

bool isInGroup(int **dataArray, int holdStop, int check)
{

	if (dataArray[3][check] - dataArray[4][holdStop] < ALLOWEDGAP)
		return true;

	if ((dataArray[1][holdStop] > dataArray[2][holdStop]) && (dataArray[1][check] > dataArray[2][check]))
	{
		if ((dataArray[2][holdStop] < dataArray[1][check]) && (dataArray[2][holdStop] > dataArray[2][check]))
			return true;
		else if ((dataArray[2][holdStop] - dataArray[1][check] > 0) && (dataArray[2][holdStop] - dataArray[1][check] < ALLOWEDMITOGAP))
			return true;

		if (dataArray[3][check] - dataArray[4][holdStop] < CONDITIONALGAP)
		{
			if (dataArray[2][holdStop] - dataArray[1][check] < CONDITIONALMITOGAP)
				return true;
		}
	}

	if ((dataArray[1][holdStop] < dataArray[2][holdStop]) && (dataArray[1][check] < dataArray[2][check]))
	{
		if ((dataArray[2][holdStop] > dataArray[1][check]) && (dataArray[2][holdStop] < dataArray[2][check]))
			return true;
		else if ((dataArray[1][check] - dataArray[2][holdStop] > 0) && (dataArray[1][check] - dataArray[2][holdStop] < ALLOWEDMITOGAP))
			return true;
		
		if (dataArray[3][check] - dataArray[4][holdStop] < CONDITIONALGAP)
		{
			if (dataArray[1][check] - dataArray[2][holdStop] < CONDITIONALMITOGAP)
				return true;
		}
	}

	return false;
}

void finalGroups(string *fileName, int sizeOfSource)
{
	ifstream input(*fileName);
	if (!input.is_open())
		return;

	bool *array = new bool[sizeOfSource];
	fill(array, array + sizeOfSource, false);

	int start, stop;

	while (!input.eof())
	{
		input >> start >> stop;
		start--;
		stop--;

		if (start < sizeOfSource && stop < sizeOfSource)
		{
			for (int i = start; i <= stop; i++)
				array[i] = true;
		}
	}
	input.close();

	ofstream output;
	int count = 0;
	bool inRegion = false;
	int holdStart, holdStop, sp = 15;
	*fileName = "final" + *fileName;
	output.open(*fileName);

	for (int i = 0; i < sizeOfSource; i++)
	{
		if (count >= ALLOWEDGAP || (i == sizeOfSource - 1 && inRegion == true))
		{
			inRegion = false;
			output << setw(sp) << holdStart + 1 << setw(sp) << holdStop + 1 << "\n";
			count = 0;
		}

		if (array[i] == false && inRegion == true)
			count++;

		if (array[i] == true && inRegion == false)
		{
			holdStart = i;
			inRegion = true;
			count = 0;
		}

		if (array[i] == true && array[i + 1] == false && inRegion == true)
		{
			holdStop = i;
			count = 0;
		}
	}

	output.close();
	delete array;
	return;
}

int regionsExtraction(string sourceFileName, string regionsFileName)
{
	ifstream sourceFile, regionsFile;
	int start = -1;
	int stop = -1;
	int pStart;
	int pStop;


	sourceFile.open(sourceFileName);
	regionsFile.open(regionsFileName);

	if (!sourceFile.is_open() || !regionsFile.is_open()){
		cout << "One or more files couldnt be opened in Regions Extraction function....";
		return 1;
	}

	bool inSeq = false;
	char hold;
	int count = 0;
	ofstream regionOut;
	string sourceName;

	//Should remove filetype from name
	if (sourceFileName.substr(sourceFileName.size() - 6, sourceFileName.size()) == ".fasta")
		sourceName = sourceFileName.substr(0, sourceFileName.size() - 6);
	else
		sourceName = sourceFileName.substr(0, sourceFileName.size() - 4);


	while (!regionsFile.eof())
	{
		pStart = start;
		pStop = stop;
		regionsFile >> start >> stop;
		if (start == pStart || stop == pStop)
			break;

		regionOut.open(sourceName + "_" + to_string(start) + " - " + to_string(stop) + ".txt");

		regionOut << ">" << sourceName + "_" + to_string(start) + " - " + to_string(stop) << "\n";

		while (count <= stop)
		{
			sourceFile.get(hold);
			if (hold == '\n' || hold == 'J')
				inSeq = true;

			if (((hold > 64 && hold < 90) || (hold > 96 && hold < 123)) && inSeq == true)
				count++;

			if (count >= start && count <= stop)
				regionOut << hold;
		}
		regionOut.close();
	}
	regionsFile.close();

	return 0;
}