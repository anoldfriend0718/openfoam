#include "fvCFD.H"
#include "keyType.H"


using namespace Foam;
int main(int argc,char* argv[]) {

    dictionary dict1;
    dict1.add(keyType("k1"),"v1");
    dict1.add(keyType("k2"),"v2");
    dict1.add(keyType("k[3-9].*",true),"v3");
    Info<<"dict1: "<<dict1<<endl;

    Info<<"dict1 toc: "<<dict1.toc()<<endl;
    Info<<"dict1 keys: "<<dict1.keys()<<endl;
    Info<<"dict1 patterns: "<<dict1.keys(true)<<endl;

    Info<<"dict1 find k1: "<<dict1.lookup("k1")<<endl;
    Info<<"dict1 find k3: "<<dict1.lookup("k3")<<endl;
    Info<<"dict1 find kk1"<<dict1.lookupOrDefault<word>("kk1","v1_new")<<endl;



    dictionary dict2;
    dict2.add(keyType("n1"),"m1");

    dictionary dict3(dict1,dict2);
    Info<<"dict3: "<<dict3<<endl;
    Info<<"dict3 parent: "<<dict3.parent()<<endl;

    dict1.add(keyType("dict2"),dict2);
    Info<<"dict1: "<<dict1<<endl;
    Info<<"dict1 subDict: "<<dict1.subDict("dict2")<<endl;

    Info<<"dict1 find kk12"<<dict1.lookupOrAddDefault<word>("kk12","v1_new")<<endl;
    Info<<"dict1: "<<dict1<<endl;
    

}