#include "fvCFD.H"


using namespace Foam;
int main(int argc,char* argv[]) {
    List<label> labelLists1{1,2};
    Info<<"labelLists1:"<<labelLists1<<" size:"<<labelLists1.size()<<nl;

    List<label> labelLists2{3,4};
    labelLists1.append(labelLists2);
    Info<<"Appended labelLists1:"<<labelLists1<<" size:"<<labelLists1.size()<<nl;

    Info<<"labelList1 max size: "<<labelLists1.max_size()<<endl;
    Info<<"labelList1 begin: "<<labelLists1.begin()<<endl;
    Info<<"labelList1 last: "<<labelLists1.last()<<endl;

    labelLists1.swap(labelLists2);
    Info<<"swap labelList1 and labelList2"<<endl;
    Info<<"Swapped labelList1: "<<labelLists1<<endl;

    reverse(labelLists1);
    Info<<"reversed labelList1: "<<labelLists1<<endl;

    sort<label>(labelLists1);
    Info<<"sorted labelList1: "<<labelLists1<<endl;

    writeListEntry(Info,labelLists1);



    return 0;
}