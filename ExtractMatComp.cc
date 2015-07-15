using namespace std;

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <dirent.h>
#include "include/ElementNames.hh"
#include "include/IsotopeMass.hh"
#include <iomanip>

// Don't forget to use GetNISTData program to generate the SetIsotopeMass() function for the IsotopeMass class
// allow GetNISTData to generate the IsotopeMass.cc file directly from a template of the class and the given NIST source file

enum  OutFilter {characters=1, numbers, NA, symbols};

void GetDataStream( string, std::stringstream&);

void FormatData(std::stringstream& stream, std::stringstream& stream2, string geoFileSourceName, bool weightPercent, string composition);
bool MovePastWord(std::stringstream& stream, string word);
string ExtractString(std::stringstream &stream, char delim, int outType=7);
void CropStream(std::stringstream& stream, int firstChar, int lastChar=0);
void FindMaterialList(std::stringstream& stream, std::vector<string> &matNameList);
void GetMaterialData(std::stringstream& stream, std::vector<string> &matNameList, std::vector<string> &matTempList, std::vector<string> &matDensList,
                    std::vector<string> &isoNameList, std::vector<double> &isoAmountVec, std::vector<double> &isoMassVec, std::vector<int> &isoStartIndex,
                    std::stringstream &original, bool wtPer);
bool FindConstructor(std::stringstream& stream, string name, std::vector<string> &matTempList, std::vector<string> &matDensList, std::vector<string> &isoNameList,
                    std::vector<double> &elemWPerVec, std::vector<double> &isoAbunVec, int &numElem, std::vector<double> &isoMassVec, std::vector<int> &elemNumIsoVec,
                    string matType, double wtPercent, double &elemMass, std::stringstream *original);
bool FindMatTempDens(std::stringstream& stream, string matName, string &matTemp, string &matDens, bool normal, std::stringstream *original);
double FindIsotopeMass(std::stringstream& stream, string isoName, std::stringstream *original);
int GetElementList(std::stringstream& stream, std::stringstream& original, std::vector<string> &matNameList, int index, std::vector<string> &isoNameList,
                std::vector<double> &elemWPerVec, std::vector<double> &isoAbunVec, int &numElem, std::vector<double> &isoMassVec, std::vector<int> &elemNumIsoVec,
                double wtPercent, std::vector<double> &wtPercentList, std::vector<double> &elemMassVec);
void FindIsotopeList(std::stringstream& stream, std::stringstream *original, string elemName, std::vector<string> &isoNameList,
                    std::vector<double> &isoAbunVec, std::vector<double> &isoMassVec, std::vector<int> &elemNumIsoVec,
                    std::vector<double> &elemWPerVec, double wtPercent, int &numElem, double &elemMass);
double GetFracComp(stringstream &stream, stringstream &original);
bool findDouble(std::stringstream *stream, string variable, double &mass);
void GetAndAddIsotope(std::stringstream& stream, std::vector<string> &isoNameList);
void GetAndAddIsotope(std::stringstream& stream, std::vector<string> &isoNameList, std::vector<double> &isoMassVec);

string CreateMacroName(string geoFileName, string outDirName);
void SetDataStream( string, std::stringstream&);



int main(int argc, char **argv)
{
    string geoFileSourceName, geoFileHeaderName, outDirName, composition;

    ElementNames elementNames;
    elementNames.SetElementNames();

    IsotopeMass isotopeMass;
    isotopeMass.SetIsotopeMass();

    std::stringstream streamS, streamH;
    string macroFileName;
    bool weightPercent=true;

    if(argc>=5&&(floor(double(argc)/2)!=ceil(double(argc)/2)))
    {
        outDirName = argv[1];
        composition = argv[2];

        if(composition=="abundance"||composition=="Abundance"||composition=="iso%"||composition=="Iso%"||composition=="isotope%"||composition=="Isotope%"||composition=="isotopic%"||composition=="Isotopic%")
        {
            weightPercent=false;
            composition="isotopic%";
        }
        else
        {
            composition="weight%";
        }
        for(int i = 3; i<argc; i+=2)
        {
            geoFileSourceName = argv[i];
            geoFileHeaderName = argv[i+1];
            GetDataStream(geoFileSourceName, streamS);
            GetDataStream(geoFileHeaderName, streamH);
            FormatData(streamS, streamH, geoFileSourceName, weightPercent, composition);
            macroFileName = CreateMacroName(geoFileSourceName, outDirName);
            SetDataStream( macroFileName, streamS);
        }
    }
    else
    {
        cout << "\nGive the the output directory, then specify weight% or abundance and then the names of the source and the header file (in that order) for each G4Stork geometry that you want to convert\n" <<  endl;
    }

    elementNames.ClearStore();
    isotopeMass.ClearStore();
}

void GetDataStream( string geoFileName, std::stringstream& ss)
{
    string* data=NULL;

    // Use regular text file
    std::ifstream thefData( geoFileName.c_str() , std::ios::in | std::ios::ate );
    if ( thefData.good() )
    {
        int file_size = thefData.tellg();
        thefData.seekg( 0 , std::ios::beg );
        char* filedata = new char[ file_size ];
        while ( thefData )
        {
            thefData.read( filedata , file_size );
        }
        thefData.close();
        data = new string ( filedata , file_size );
        delete [] filedata;
    }
    else
    {
    // found no data file
    //                 set error bit to the stream
        ss.setstate( std::ios::badbit );
    }
    if (data != NULL)
    {
        ss.str(*data);
        if(data->back()!='\n')
            ss << "\n";
        ss.seekg( 0 , std::ios::beg );
    }

    delete data;
}

void FormatData(std::stringstream& stream, std::stringstream& stream2, string geoFileSourceName, bool weightPercent, string composition)
{
    stringstream numConv;

    std::vector<string> matNameList;
    matNameList.reserve(20);
    std::vector<string> matTempList;
    matNameList.reserve(20);
    std::vector<string> matDensList;
    matNameList.reserve(20);
    std::vector<string> isoNameList;
    isoNameList.reserve(200);
    std::vector<double> isoAmountVec;
    isoAmountVec.reserve(200);
    std::vector<double> isoMassVec;
    isoMassVec.reserve(200);
    std::vector<int> isoStartIndex;
    isoStartIndex.reserve(20);
    std::stringstream original;

    original.str(stream2.str()+stream.str());

    MovePastWord(stream, "::ConstructMaterials()");
    int pos = stream.tellg();
    CropStream(stream, pos);
    FindMaterialList(stream, matNameList);

    GetMaterialData(stream, matNameList, matTempList, matDensList, isoNameList, isoAmountVec, isoMassVec, isoStartIndex, original, weightPercent);

    stream.str("");
    stream.clear();

    stream << "This file lists the " << composition << " and the mass in amu, of all the isotopes used to create each material used in " << geoFileSourceName << endl;

    isoStartIndex.push_back(isoNameList.size());

    for(int i=0; i<int(matNameList.size()); i++)
    {
        stream.fill('-');
        stream << std::setw(84) << std::left << "Material: "+matNameList[i] << '\n' << endl;
        stream << "Density: " << matDensList[i] << '\n' << "Number of Isotopes: " << isoStartIndex[i+1]-isoStartIndex[i] << '\n' << endl;

        stream.fill(' ');
        stream << std::setw(20) << std::left << "Isotope Name:" << std::setw(20) << std::left << "Amount:" << std::setw(20) << std::left << "Mass:"
                  << std::setw(20) << std::left << "Temperature:" << '\n' << endl;

        for(int j=isoStartIndex[i]; j<isoStartIndex[i+1]; j++)
        {
            numConv.str("");
            numConv.clear();
            numConv << isoAmountVec[j];
            stream << std::setw(20) << std::left << isoNameList[j] << std::setw(20) << std::left << numConv.str()+composition;

            numConv.str("");
            numConv.clear();
            numConv << isoMassVec[j];
            stream << std::setw(20) << std::left << numConv.str()+"amu" << matTempList[i] << '\n';
        }

        stream << "\n";
    }
}

bool MovePastWord(std::stringstream& stream, string word)
{
    std::vector<string> wordParts;
    int pos=0, start;
    bool check=true, firstPass=true;

    start = stream.tellg();

    for(int i=0; i<int(word.length()); i++)
    {
        if(word[i]==' ')
        {
            if(check)
            {
                pos=i+1;
            }
            else
            {
                wordParts.push_back(word.substr(pos,i-pos));
                pos=i+1;
                check=true;
            }
        }
        else
        {
            check=false;
            if(i==int(word.length()-1))
            {
                wordParts.push_back(word.substr(pos,i-pos+1));
            }
        }
    }

    if(wordParts.size()==0)
    {
        wordParts.push_back(word);
    }

    string wholeWord, partWord;
    check=false;
    char line[256];

    while(!check)
    {
        if(!stream)
        {
            if(firstPass)
            {
                stream.clear();
                stream.seekg(start, std::ios::beg);
                firstPass=false;
            }
            else
            {
                break;
            }
        }
        if(stream.peek()=='/')
        {
            stream.get();
            if(stream.peek()=='/')
            {
                stream.getline(line,256);
            }
            else if(stream.peek()=='*')
            {
                stream.get();
                while(stream)
                {
                    if(stream.get()=='*')
                    {
                        if(stream.get()=='/')
                        {
                            break;
                        }
                    }
                }
            }
        }
        else if(stream.peek()=='\n')
        {
            stream.getline(line,256);
        }
        else if(stream.peek()=='\t')
        {
            stream.get();
        }
        else if(stream.peek()==' ')
        {
            stream.get();
        }
        else
        {
            for(int i=0; i<int(wordParts.size()); i++)
            {
                stream >> wholeWord;
                if(int(wholeWord.length())>=int((wordParts[i]).length()))
                {
                    if(firstPass)
                    {
                        check=(wholeWord==(wordParts[i]));
                        if(!check)
                        {
                            break;
                        }
                    }
                    else
                    {
                        partWord = wholeWord.substr(0, (wordParts[i]).length());
                        check=(partWord==(wordParts[i]));

                        if(check)
                        {
                            stream.seekg((partWord.length()-wholeWord.length()),std::ios_base::cur);
                        }
                        else if(0==i)
                        {
                            partWord = wholeWord.substr(wholeWord.length()-(wordParts[i]).length(), (wordParts[i]).length());
                            check=(partWord==(wordParts[i]));
                        }

                        if(!check)
                        {
                            break;
                        }
                    }

                }
                else
                {
                    break;
                }
            }
        }

    }

    if(!check)
    {
        stream.clear();
        stream.seekg(start, std::ios::beg);
    }

    return check;
}

string ExtractString(std::stringstream &stream, char delim, int outType)
{
    string value="";
    bool charOut=false, numOut=false, symOut=false;
    char letter;
    //bool first=true;

    if(outType==0)
    {

    }
    else if(outType==1)
    {
        charOut=true;
    }
    else if(outType==2)
    {
        numOut=true;
    }
    else if(outType==3)
    {
        charOut=true;
        numOut=true;
    }
    else if(outType==4)
    {
        symOut=true;
    }
    else if(outType==5)
    {
        charOut=true;
        symOut=true;
    }
    else if(outType==6)
    {
        numOut=true;
        symOut=true;
    }
    else
    {
        charOut=true;
        numOut=true;
        symOut=true;
    }

    while(stream&&(stream.peek()!=delim))
    {
        letter = stream.get();
        if(((letter>='A')&&(letter<='Z'))||((letter>='a')&&(letter<='z')))
        {
            if(charOut)
            {
                value+=letter;
                //first=true;
            }
            /*else if(first)
            {
                value+=' ';
                first=false;
            }*/
        }
        else if(((letter>='0')&&(letter<='9'))||(letter=='.')||(letter=='-'))
        {
            if(numOut)
            {
                value+=letter;
                //first=true;
            }
            /*else if(first)
            {
                value+=' ';
                first=false;
            }*/
        }
        else
        {
            if(symOut)
            {
                value+=letter;
                //first=true;
            }
            /*else if(first)
            {
                value+=' ';
                first=false;
            }*/
        }
    }
    return value;
}

void CropStream(std::stringstream& stream, int firstChar, int lastChar)
{
    if(lastChar==0)
        stream.seekg( 0 , std::ios::end );
    else
        stream.seekg( lastChar , std::ios::beg );

    int file_size = int(stream.tellg())-firstChar;
    stream.seekg( firstChar , std::ios::beg );
    char* filedata = new char[ file_size ];

    stream.read( filedata , file_size );
    if(!stream)
    {
        cout << "\n #### Error reading string stream ###" << endl;
        return;
    }
    stream.str("");
    stream.clear();

    stream.write( filedata, file_size);
    stream.seekg(0 , std::ios::beg);

    delete filedata;
}

void FindMaterialList(std::stringstream& stream, std::vector<string> &matNameList)
{
    string name="";
    while(MovePastWord(stream, "matMap["))
    {
        stream.get();

        ExtractString(stream, '=', 0);

        stream.get();

        name=ExtractString(stream, ';', int(characters+numbers));
        if(name!="")
        {
            matNameList.push_back(name);
        }
        else
        {
            cout << "\nError: found a blank when trying to extract material name\n" << endl;
        }
        name.clear();

    }
    stream.clear();
    stream.seekg(0, std::ios::beg);
}

void GetMaterialData(std::stringstream& stream, std::vector<string> &matNameList, std::vector<string> &matTempList, std::vector<string> &matDensList,
                    std::vector<string> &isoNameList, std::vector<double> &isoAmountVec, std::vector<double> &isoMassVec, std::vector<int> &isoStartIndex,
                    std::stringstream &original, bool wtPer)
{
    std::vector<double> elemWPerVec;
    std::vector<double> isoAbunVec;
    std::vector<double> wtPercentList;
    std::vector<double> elemMassVec;
    std::vector<double> wtMassFracVec;
    std::vector<int> elemVecSizeList;
    elemVecSizeList.reserve(100);
    std::vector<int> elemNumIsoVec;
    elemNumIsoVec.reserve(100);

    double wtPercent, elemMass, sumWMassFrac=0, sumAmount;
    int addMat=0, extraMat=0, totalAddMat=0, numElem, isoIndex=0, isoIndex2, elemIndex=0;

    for(int i=0; i<int(matNameList.size()); i++)
    {
        if(addMat==0)
        {
            wtPercent=1.0;
            numElem=0;
            sumWMassFrac=0.;
            totalAddMat=0;
            elemMassVec.clear();
            wtPercentList.clear();
            elemWPerVec.clear();
            isoAbunVec.clear();
            wtMassFracVec.clear();
        }
        else
        {
            addMat--;
            wtPercent=wtPercentList[addMat];
        }

        elemMass=0.;
        stream.clear();
        stream.seekg(0, std::ios::beg);
        if(FindConstructor(stream, matNameList[i], matTempList, matDensList, isoNameList, elemWPerVec, isoAbunVec, numElem, isoMassVec, elemNumIsoVec, "Material", wtPercent, elemMass, &original))
        {
            extraMat=GetElementList(stream, original, matNameList, i, isoNameList, elemWPerVec, isoAbunVec, numElem, isoMassVec, elemNumIsoVec, wtPercent, wtPercentList, elemMassVec);
            addMat+=extraMat;
            totalAddMat+=extraMat;
        }
        else
        {
            if(numElem>int(elemMassVec.size()))
                elemMassVec.push_back(elemMass);
        }
        stream.clear();
        stream.seekg(0, std::ios::beg);

        if(addMat==0)
        {
            elemVecSizeList.push_back(numElem);
            matNameList.erase(matNameList.begin()+i+1-totalAddMat, matNameList.begin()+i+1);
            matTempList.erase(matTempList.begin()+i+1-totalAddMat, matTempList.begin()+i+1);
            matDensList.erase(matDensList.begin()+i+1-totalAddMat, matDensList.begin()+i+1);
            i=i-totalAddMat;
            sumAmount=0;

            if(wtPer)
            {
                isoIndex2=isoIndex;
                isoStartIndex.push_back(isoIndex2);
                for(int k=elemIndex; k<int(elemNumIsoVec.size()); k++)
                {
                    for(int j=isoIndex; j<(isoIndex+elemNumIsoVec[k]); j++)
                    {
                        isoAmountVec.push_back(isoMassVec[j]*isoAbunVec[j-isoIndex2]*elemWPerVec[k-elemIndex]/elemMassVec[k-elemIndex]);
                        sumAmount+=isoAmountVec.back();
                    }
                    isoIndex=isoIndex+elemNumIsoVec[k];
                }
                elemIndex=elemNumIsoVec.size();
            }
            else
            {

                for(int k=0; k<int(elemVecSizeList.back()); k++)
                {
                    wtMassFracVec.push_back(elemWPerVec[k]/elemMassVec[k]);
                    sumWMassFrac+=wtMassFracVec[k];
                }

                isoIndex2=isoIndex;
                isoStartIndex.push_back(isoIndex2);
                for(int k=0; k<int(elemVecSizeList.back()); k++)
                {
                    wtMassFracVec[k]=wtMassFracVec[k]/sumWMassFrac;

                    for(int j=isoIndex; j<(isoIndex+elemNumIsoVec[k+elemIndex]); j++)
                    {
                        isoAmountVec.push_back(isoAbunVec[j-isoIndex2]*wtMassFracVec[k]);
                        sumAmount+=isoAmountVec.back();
                    }
                    isoIndex=isoIndex+elemNumIsoVec[k+elemIndex];
                }
                elemIndex=elemNumIsoVec.size();
            }

            for(int k=isoIndex2; k<int(isoNameList.size()); k++)
            {
                for(int j=k+1; j<int(isoNameList.size()); j++)
                {
                    if(isoNameList[j]==isoNameList[k])
                    {
                        isoAmountVec[k]+=isoAmountVec[j];
                        isoNameList.erase(isoNameList.begin()+j);
                        isoAmountVec.erase(isoAmountVec.begin()+j);
                        isoMassVec.erase(isoMassVec.begin()+j);
                        isoIndex--;
                        j--;
                    }
                }
            }

            for(int j=isoIndex2; j<int(isoAmountVec.size()); j++)
            {
                isoAmountVec[j]/=sumAmount;
            }
        }
    }
}

bool FindConstructor(std::stringstream& stream, string name, std::vector<string> &matTempList, std::vector<string> &matDensList, std::vector<string> &isoNameList,
                    std::vector<double> &elemWPerVec, std::vector<double> &isoAbunVec, int &numElem, std::vector<double> &isoMassVec, std::vector<int> &elemNumIsoVec,
                    string matType, double wtPercent, double &elemMass, std::stringstream *original)
{
    string check, matTemp, matDens;
    std::stringstream line;
    int pos;
    if(!MovePastWord(stream, (name+" =")))
    {
        cout << "\nError: could not find constructor for " << name << "\n" << endl;
        return false;
    }
    else
    {
        pos=stream.tellg();
        int count=0;

        if(matType=="Material")
        {
            bool intType=true;
            int pos1=pos;
            while(stream.peek()!=')')
            {
                if(stream.get()==',')
                {
                    if(count==0)
                    {
                        pos1=stream.tellg();
                    }
                    count++;
                }
                if(count==3)
                {
                    line.str(ExtractString(stream, ')', int(characters+numbers+symbols))+")");
                    line.seekg(0,std::ios::beg);
                    bool test=false;
                    while((line.peek()!=')')&&(line.peek()!=','))
                    {
                        if(line.get()=='.')
                        {
                            intType=false;
                            break;
                        }
                        if(line.peek()==',')
                        {
                            test=true;
                        }
                    }
                    if(intType)
                    {
                        line.seekg(0,std::ios::beg);
                        string variable;
                        double num;
                        if(test)
                        {
                            variable = ExtractString(line, ',', int(characters+numbers));
                        }
                        else
                        {
                            variable = ExtractString(line, ')', int(characters+numbers));
                        }
                        if (findDouble(original, variable, num))
                        {
                            stringstream numConv;
                            numConv << num;
                            while(numConv)
                            {
                                if(numConv.get()=='.')
                                {
                                    intType=false;
                                    break;
                                }
                            }
                        }
                    }
                    if(intType)
                    {
                        intType = (!MovePastWord(line, "kState"));
                    }
                    line.str("");
                    line.clear();
                    break;
                }
            }
            stream.seekg(pos1, std::ios::beg);
            if(intType)
            {
                FindMatTempDens(stream, name, matTemp, matDens, true, original);
                matTempList.push_back(matTemp);
                matDensList.push_back(matDens);
                stream.seekg(pos1, std::ios::beg);
                return true;
            }
            else
            {
                FindMatTempDens(stream, name, matTemp, matDens, false, original);
                matTempList.push_back(matTemp);
                matDensList.push_back(matDens);
                stream.seekg(pos1, std::ios::beg);
                elemWPerVec.push_back(wtPercent);
                isoAbunVec.push_back(1.0);
                elemNumIsoVec.push_back(1);
                GetAndAddIsotope(stream, isoNameList, isoMassVec);
                elemMass=isoMassVec.back();
                numElem++;
                return false;
            }
        }
        else if(matType=="Element")
        {
            int count=0, pos1=0, pos2=0;
            while(stream.peek()!=';')
            {
                if(stream.get()==',')
                {
                    pos1=pos2;
                    pos2=stream.tellg();
                    count++;
                }
            }
            if(count==3)
            {
                stream.seekg(pos1, std::ios::beg);
                isoAbunVec.push_back(1.0);
                elemNumIsoVec.push_back(1);
                GetAndAddIsotope(stream, isoNameList, isoMassVec);
                elemMass=isoMassVec.back();
                return false;
            }
            else
            {
                stream.seekg(pos, std::ios::beg);
                return true;
            }
        }
        else
        {
            return true;
        }
    }
}

bool FindMatTempDens(std::stringstream& stream, string matName, string &matTemp, string &matDens, bool normal, std::stringstream *original)
{
    int count=0, limit, limit2;
    bool standard=true, number=true, first=true, found=false;
    matTemp="";
    matDens="";
    double temperature;
    char letter;
    std::stringstream temp;

    if(normal)
    {
        limit=3;
        limit2=0;
    }
    else
    {
        limit=4;
        limit2=2;
    }

    while(true)
    {
        letter = stream.get();
        if(letter==',')
        {
            count++;
        }
        if(limit2==count)
        {
            break;
        }
    }

    while((stream.peek()!=',')&&(stream.peek()!=')'))
    {
        letter=stream.get();
        if(((letter>='0')&&(letter<='9'))||(letter=='.')||(letter=='-'))
        {
            if(first)
            {
                number=true;
                first=false;
            }
            matDens += letter;
        }
        else if((((letter>='A')&&(letter<='Z'))||((letter>='a')&&(letter<='z')))||(letter=='[')||(letter==']')||(letter==',')||(letter=='/'))
        {
            if(first)
            {
                number=false;
                first=false;
            }
            matDens += letter;
        }
    }

    if(matDens!="")
    {
        if(number)
        {
            found=true;
        }
        else if(original!=NULL)
        {
            found=findDouble(original, matDens, temperature);
            temp << temperature;
            temp >> matDens;
            matDens = matDens+"g/cm3";
            temp.str("");
            temp.clear();
        }
    }
    else
    {
        cout << "\nError: unable to find density for " << matName << " in the expected position\n" << endl;
    }

    number=first=true;

    while(stream.peek()!=';')
    {
        if((stream.get())==',')
        {
            count++;
        }
        if(limit==count)
        {
            standard=false;
            break;
        }
    }
    if(standard)
    {
        matTemp="273.15kelvin";
    }
    else
    {
        while((stream.peek()!=',')&&(stream.peek()!=')'))
        {
            letter=stream.get();
            if(((letter>='0')&&(letter<='9'))||(letter=='.')||(letter=='-'))
            {
                if(first)
                {
                    number=true;
                    first=false;
                }
                matTemp+= letter;
            }
            else if((((letter>='A')&&(letter<='Z'))||((letter>='a')&&(letter<='z')))||(letter=='[')||(letter==']')||(letter==',')||(letter=='/'))
            {
                if(first)
                {
                    number=false;
                    first=false;
                }
                matTemp+= letter;
            }
        }

        if(matTemp!="")
        {
            if(number)
            {
                found=true;
            }
            else if(original!=NULL)
            {
                found=findDouble(original, matTemp, temperature);
                temp << temperature;
                temp >> matTemp;
                matTemp = matTemp+"kelvin";
                temp.str("");
                temp.clear();
            }
        }
        else
        {
            cout << "\nError: unable to find temperature for " << matName << " in the expected position\n" << endl;
        }
    }

    return found;
}

int GetElementList(std::stringstream& stream, std::stringstream& original, std::vector<string> &matNameList, int index, std::vector<string> &isoNameList,
                std::vector<double> &elemWPerVec, std::vector<double> &isoAbunVec, int &numElem, std::vector<double> &isoMassVec, std::vector<int> &elemNumIsoVec,
                double wtPercent, std::vector<double> &wtPercentList, std::vector<double> &elemMassVec)
{
    string name="";
    int addMat=0, count=0, pos;
    std::stringstream checkCon, numConv;
    double tempNum, elemMass, fracSum=0;
    bool abun=true;

    while(MovePastWord(stream, matNameList[index]+" ->"))
    {
        elemMass=0.;

        name.clear();
        name=ExtractString(stream, '(', int(characters));

        stream.get();

        if(name=="AddElement")
        {
            name.clear();
            name=ExtractString(stream, ',', int(characters+numbers));
            stream.get();
            checkCon.clear();
            checkCon.str(name);
            if(MovePastWord(checkCon, "new G4Element"))
            {
                ExtractString(stream, ',', 0);
                stream.get();
                GetAndAddIsotope(stream, isoNameList, isoMassVec);
                elemMass=isoMassVec.back();
                elemMassVec.push_back(elemMass);
                ExtractString(stream, ')', 0);
                ExtractString(stream, ',', 0);
                stream.get();

                pos=stream.tellg();
                numConv.str(ExtractString(stream, ')', int(numbers)));

                while(numConv)
                {
                    if(numConv.get()=='.')
                    {
                        abun=false;
                    }
                }

                stream.seekg(pos, std::ostream::beg);
                tempNum=GetFracComp(stream,original);
                fracSum+=tempNum*elemMass;

                numConv.clear();
                numConv.seekg(0, std::ostream::beg);

                elemWPerVec.push_back(wtPercent*tempNum);
                isoAbunVec.push_back(1.0);
                elemNumIsoVec.push_back(1);
                numElem++;
                count++;
            }
            else if(name!="")
            {
                pos=stream.tellg();
                numConv.str(ExtractString(stream, ')', int(numbers)));

                while(numConv)
                {
                    if(numConv.get()=='.')
                    {
                        abun=false;
                    }
                }
                stream.seekg(pos, std::ostream::beg);
                tempNum=GetFracComp(stream,original);

                numConv.clear();
                numConv.seekg(0, std::ostream::beg);

                elemWPerVec.push_back(wtPercent*tempNum);
                numElem++;
                count++;

                FindIsotopeList(stream, &original, name, isoNameList, isoAbunVec, isoMassVec, elemNumIsoVec, elemWPerVec, elemWPerVec.back(), numElem, elemMass);
                elemMassVec.push_back(elemMass);
                fracSum+=tempNum*elemMass;

            }
            else
            {
                cout << "\nError: found a blank when trying to extract element name\n" << endl;
            }
        }
        else if(name=="AddMaterial")
        {
            name.clear();
            abun=false;

            name=ExtractString(stream, ',', int(characters+numbers));
            stream.get();
            checkCon.clear();
            checkCon.str(name);
            if(MovePastWord(checkCon, "new G4Material"))
            {
                GetAndAddIsotope(stream, isoNameList, isoMassVec);
                ExtractString(stream, ')', 0);
                ExtractString(stream, ',', 0);
                stream.get();
                tempNum=GetFracComp(stream,original);
                elemWPerVec.push_back(wtPercent*tempNum);
                isoAbunVec.push_back(1.0);
                elemNumIsoVec.push_back(1);
                numElem++;
                count++;
            }
            else if(name!="")
            {
                matNameList.insert(matNameList.begin()+index+1,name);
                addMat++;
                tempNum=GetFracComp(stream,original);
                wtPercentList.push_back(wtPercent*tempNum);
            }
            else
            {
                cout << "\nError: found a blank when trying to extract material name\n" << endl;
            }
        }
    }

    if(abun)
    {
        for(int i=numElem-count; i<numElem; i++)
        {
            elemWPerVec[i]*=elemMassVec[i]/fracSum;
        }
    }

    return addMat;
}

void FindIsotopeList(std::stringstream& stream, std::stringstream *original, string elemName, std::vector<string> &isoNameList, std::vector<double> &isoAbunVec,
                    std::vector<double> &isoMassVec, std::vector<int> &elemNumIsoVec, std::vector<double> &elemWPerVec, double wtPercent, int &numElem, double &elemMass)
{
    std::stringstream checkCon, numConv;
    int numIso=0, pos, pos2, isoIndex=isoNameList.size(), isoIndex2=isoAbunVec.size();
    double tempNum;
    std::vector<string> ignore;

    pos2=stream.tellg();

    stream.clear();
    stream.seekg(0, std::ios::beg);

    if(FindConstructor(stream, elemName, ignore, ignore, isoNameList, elemWPerVec, isoAbunVec, numIso, isoMassVec, elemNumIsoVec, "Element", wtPercent, elemMass, original))
    {
        string name="";
        while(MovePastWord(stream, elemName+" ->"))
        {
            name.clear();
            name=ExtractString(stream, '(', int(characters));

            stream.get();

            if(name=="AddIsotope")
            {
                name.clear();
                name=ExtractString(stream, ',', int(characters+numbers));
                stream.get();
                checkCon.clear();
                checkCon.str(name);
                if(MovePastWord(checkCon, "new G4Isotope"))
                {
                    string test, test2;
                    GetAndAddIsotope(stream, isoNameList);
                    isoMassVec.push_back(FindIsotopeMass(stream, isoNameList.back(), original));
                    ExtractString(stream, ')', 0);
                    ExtractString(stream, ',', 0);
                    stream.get();

                    tempNum=GetFracComp(stream, *original);
                    isoAbunVec.push_back(tempNum);
                    numIso++;
                }
                else if(name!="")
                {
                    tempNum=GetFracComp(stream, *original);
                    pos=stream.tellg();
                    stream.clear();
                    stream.seekg(0, std::ios::beg);

                    if(FindConstructor(stream, name, ignore, ignore, isoNameList, elemWPerVec, isoAbunVec, numIso, isoMassVec, elemNumIsoVec, "Isotope", wtPercent, elemMass, original))
                    {
                        ExtractString(stream, ',', 0);
                        stream.get();
                        GetAndAddIsotope(stream, isoNameList);
                        isoMassVec.push_back(FindIsotopeMass(stream, isoNameList.back(), original));

                        isoAbunVec.push_back(tempNum);
                        numIso++;
                        stream.seekg(pos, std::ios::beg);
                    }
                    else
                    {
                        stream.seekg(pos, std::ios::beg);
                    }
                }
                else
                {
                    cout << "\nError: found a blank when trying to extract isotope name\n" << endl;
                }
            }
        }
    }

    elemNumIsoVec.push_back(numIso);
    double sum=0;

    for(int i=0; i<elemNumIsoVec.back(); i++)
    {
        elemMass+=(isoMassVec[isoIndex+i])*(isoAbunVec[isoIndex2+i]);
        sum+=isoAbunVec[isoIndex2+i];
    }

    elemMass/=sum;

    for(int i=0; i<elemNumIsoVec.back(); i++)
    {
        isoAbunVec[isoIndex2+i]=isoAbunVec[isoIndex2+i]/sum;
    }

    stream.seekg(pos2, std::ios::beg);

}

double GetFracComp(stringstream &stream, stringstream &original)
{
    string test;
    int pos, endPos=0;
    bool number=true;
    char letter;
    stringstream numConv;
    double tempNum;

    pos=stream.tellg();
    numConv.str(ExtractString(stream, ')', int(numbers)));
    numConv >> tempNum;
    stream.seekg(pos, std::ios::beg);
    test=ExtractString(stream, ')', int(numbers+symbols+characters));
    numConv.str(test);

    for(int i=0; i<int(test.length()); i++)
    {
        letter = test[i];
        if(((letter>='A')&&(letter<='Z'))||((letter>='a')&&(letter<='z')))
        {
            number=false;
            break;
        }
        else if((letter>='1')&&(letter<='9'))
        {
            number=true;
            break;
        }
        else if((letter=='+')||(letter=='-')||(letter=='*')||(letter=='/'))
        {
            endPos=i;
        }

    }

    if(endPos!=0)
    {
        test=test.substr(0,endPos+1);
    }

    if(!number)
    {
        findDouble(&original, test, tempNum);
    }

    numConv.clear();
    numConv.seekg(0, std::iostream::beg);

    if(MovePastWord(numConv, "perCent"))
        tempNum=tempNum/100;

    return tempNum;

}

double FindIsotopeMass(std::stringstream& stream, string isoName, std::stringstream *original)
{
    bool number=true, first=true;
    double mass=0.;
    char letter;
    std::stringstream temp;

    while((stream.peek()!=',')&&(stream.peek()!=')'))
    {
        letter=stream.get();
        if(((letter>='0')&&(letter<='9'))||(letter=='.')||(letter=='-'))
        {
            if(first)
            {
                number=true;
                first=false;
            }
            temp << letter;
        }
        else if((((letter>='A')&&(letter<='Z'))||((letter>='a')&&(letter<='z')))||(letter=='[')||(letter==']')||(letter==','))
        {
            if(first)
            {
                number=false;
                first=false;
            }
            if(!number)
                temp << letter;
        }
    }

    if(temp.str()!="")
    {
        if(number)
        {
            temp >> mass;
        }
        else if(original!=NULL)
        {
            findDouble(original, temp.str(), mass);
        }
    }
    else
    {
        cout << "\nError: unable to find mass for " << isoName << " in the expected position\n" << endl;
    }

    return mass;
}

bool findDouble(std::stringstream *stream, string variable, double &mass)
{
    bool arrayElem=false, number=false, first=true;
    std::vector<int> arrayIndex;
    int index, pos1, pos2, count=0;
    stringstream numConv, temp;
    char letter;
    stream->seekg(0, std::ios::beg);

    while(variable.back()==']')
    {
        arrayElem=true;
        pos1=variable.find_first_of('[',0);
        pos2=variable.find_first_of(']',0);
        numConv.str(variable.substr(pos1,pos1-pos2-1));
        numConv >> index;
        numConv.clear();
        numConv.str("");
        arrayIndex.push_back(index);
        variable.erase(pos1, pos2-pos1+1);
    }

    if(variable!="")
    {
        if(*stream)
        {
            if(MovePastWord((*stream),variable+" ="))
            {
                temp.str(ExtractString((*stream),';',int(numbers+characters)));
                temp.str(temp.str()+';');

                if(arrayElem)
                {
                    for(int i=0; i<int(arrayIndex.size()); i++)
                    {
                        ExtractString(temp,'{',0);
                        temp.get();
                        while(count!=arrayIndex[i])
                        {
                            letter=temp.get();
                            if(letter=='{')
                            {
                               ExtractString(temp,'}',0);
                               temp.get();
                            }
                            else if(letter==',')
                            {
                                count++;
                            }
                        }
                    }
                }
                while((temp.peek()!=',')&&(temp.peek()!=';'))
                {
                    letter=temp.get();
                    if(((letter>='0')&&(letter<='9'))||(letter=='.')||(letter=='-'))
                    {
                        if(first)
                        {
                            number=true;
                            first=false;
                        }
                        numConv << letter;
                    }
                    else if((((letter>='A')&&(letter<='Z'))||((letter>='a')&&(letter<='z')))||(letter=='[')||(letter==']')||(letter==','))
                    {
                        if(first)
                        {
                            number=false;
                            first=false;
                        }
                        if(!number)
                            numConv << letter;
                    }
                }
                if(number)
                {
                    numConv >> mass;
                }
                else
                {
                    (*stream).seekg(0,std::ios::beg);
                    return findDouble(stream, numConv.str(), mass);
                }

            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        return false;
    }

    return true;
}

void GetAndAddIsotope(std::stringstream& stream, std::vector<string> &isoNameList)
{
    std::stringstream isoName;
    ElementNames* elementNames;
    int Z;

    isoName << ExtractString(stream, ',', int(numbers));

    stream.get();
    isoName >> Z;
    isoName.clear();
    isoName << '_';

    isoName << (ExtractString(stream, ',', int(numbers))).c_str();
    stream.get();
    isoName << "_" << elementNames->GetName(Z);
    isoNameList.push_back(isoName.str());

    isoName.str("");
    isoName.clear();
}

void GetAndAddIsotope(std::stringstream& stream, std::vector<string> &isoNameList, std::vector<double> &isoMassVec)
{
    std::stringstream isoName, numConv;
    ElementNames* elementNames;

    IsotopeMass *isotopeMass;

    int Z, A;
    double mass;

    isoName << ExtractString(stream, ',', int(numbers));

    stream.get();
    isoName >> Z;
    isoName.clear();
    isoName << '_';

    numConv << (ExtractString(stream, ',', int(numbers))).c_str();
    stream.get();
    isoName << (numConv.str()).c_str();
    numConv >> A;
    isoName << "_" << elementNames->GetName(Z);
    mass = isotopeMass->GetIsotopeMass(Z,A);
    isoNameList.push_back(isoName.str());
    isoMassVec.push_back(mass);

    isoName.str("");
    isoName.clear();
}

string CreateMacroName(string geoFileName, string outDirName)
{
    if((geoFileName.substr(geoFileName.length()-3,3))==".cc")
    {
        geoFileName=geoFileName.substr(0,geoFileName.length()-3);
    }
    size_t pos = geoFileName.find_last_of('/');
    size_t pos2 = std::string::npos;
    if(pos == std::string::npos)
        pos=0;
    else
        pos++;

    if(geoFileName.length()>11)
    {
        string test = geoFileName.substr(geoFileName.length()-11, 11);
        if((test=="Constructor")||(test=="constructor"))
        {
            pos2 = geoFileName.length()-11;
        }
    }

    return (outDirName+"MatComp"+geoFileName.substr(pos, pos2-pos)+".txt");
}

void SetDataStream( string macroFileName, std::stringstream& ss)
{
  std::ofstream out( macroFileName.c_str() , std::ios::out | std::ios::trunc );
  if ( ss.good() )
  {
     ss.seekg( 0 , std::ios::end );
     int file_size = ss.tellg();
     ss.seekg( 0 , std::ios::beg );
     char* filedata = new char[ file_size ];

     while ( ss )
     {
        ss.read( filedata , file_size );
        if(!file_size)
        {
            cout << "\n #### Error the size of the stringstream is invalid ###" << endl;
            break;
        }
     }

     out.write(filedata, file_size);
     if (out.fail())
    {
        cout << endl << "writing the ascii data to the output file " << macroFileName << " failed" << endl
             << " may not have permission to delete an older version of the file" << endl;
    }
     out.close();
     delete [] filedata;
  }
  else
  {
// found no data file
//                 set error bit to the stream
     ss.setstate( std::ios::badbit );

     cout << endl << "### failed to write to ascii file " << macroFileName << " ###" << endl;
  }
   ss.str("");
}
