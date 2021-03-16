#include "fvCFD.H"
#include "messageStream.H"
#include "stringOps.H"


using namespace Foam;
int main(int argc,char* argv[]) {

    //string:  A class for handling character strings derived from std::string.

    // Strings may contain any characters and therefore are delimited by quotes
    // for IO : "any list of characters".

    // Used as a base class for word and fileName.

    string s1("~OpenFOAM/controlDict");

    Info<<"s1 expand: "<<s1.expand()<<endl;
    Info<<"s1 hash: "<<string::hash()(s1)<<endl;

    s1.replace("controlDict","controlDict1");
    Info<<"s1 replaced: "<<s1<<endl;

    string::size_type i=s1.find("control");
    if(i!=string::npos)
    {
        Info<<"found control: "<<i<<endl;
    }

    string s2="   s2";
    Info<<"s2 trimmed: "<<stringOps::trimLeft(s2)<<endl;

    string s3="ssdd";
    Info<<"s3 isremoveRepeated: "<<s3.removeRepeated('d')<<" s3: "<<s3<<endl;


    IOobject::writeDivider(Info);

    word::debug=1;
    word w1("{'/w  e}");
    Info<<"w1: "<<w1<<endl; //when debug=1,strip invalid character, like " \ / ; {} and whiteSpace



    return 0;
}