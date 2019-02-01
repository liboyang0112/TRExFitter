// The scripts behave exactly as "hadd -f", with the notable difference 
// that it will not add identical TH1, keeping just one copy of them
// So if calling
//     hupdate targetfile file1 file2
// with
//     file1 containing two TH1:   h1 h2
//     file2 containing two TH1:   h2 h3
// targetfile will contain:        h1 h2 h3, 
// with h2 taken only from file1 (not added to its copy in file2)
// The script is not (yet) iterating over subdirectories, not needed for its purposes
//
// author: Giancarlo Panizzo giancarlo.panizzo@cern.ch
//         (adapted from hadd)
#include <stdio.h>
#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

void MergeRootfile( TDirectory *target, TList *sourcelist );

// main function
// -------------------------------------------------------
// -------------------------------------------------------
TList *FileList;
TFile *Target;
TFile*whichfile;

int main(int argc, char **argv){
   
   std::cout << "Creating target mergin histograms from source files. Histograms not summed if have the same name: only the first one is kept." << std::endl;
   
   std::cout << "Target file:\t" << argv[1] << std::endl;
   
   Target = new TFile( argv[1], "RECREATE" );
   FileList = new TList();
   
   TDirectory *dir = gDirectory;
      
   if (argc<3) {
      std::cout << "Too few arguments ..." <<std::endl    ;
      return -1;
   }
  
   for (int i=2; i<argc; i++){
      std::cout << "Source file("<< i-1<<"):\t" << argv[i] <<std::endl;
      whichfile=new TFile(argv[i]);
      FileList->Add( whichfile );
      dir->cd();
   }

  // call the function
  MergeRootfile( Target, FileList );

  Target->Close();
  
  return 0;
}


void MergeRootfile( TDirectory *target, TList *sourcelist ) {
   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove( 0, 2 );

   TList *ListAllKeys=new TList();  // This is supposed to keep track of all different objects in all files...
 
   // take next file:
   TIter nextsource_iter( sourcelist );
//    TDirectory *current_sourcedir = gDirectory;

   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   TFile *next_source; 
   while ( (next_source = (TFile*)nextsource_iter() ) ) {//  start cycle over files
      next_source->cd( path );
      // loop over all keys in this directory (meaning of the directory of the first file! This is the point in which this script starts to be different from hadd)
      TChain *globChain = 0;
      TIter nextkey( next_source->GetListOfKeys() );
      TKey *key, *oldkey=0;
      while ( (key = (TKey*)nextkey())) {
        //keep only the highest cycle number for each key
         if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue; 
	   // ok, this key is the only one in -this- file. 
           // Then check if it was already present in other files.
   	   TIter nextkeyAllFiles(ListAllKeys);
	   Bool_t skipkey=kFALSE;
	   TKey* keyGlobal;
           while ( (keyGlobal = (TKey*)nextkeyAllFiles()) ) if ( !strcmp(keyGlobal->GetName(),key->GetName())) skipkey=kTRUE; //FIXME: no check if the histograms are actually different, just look at the name ...
	   if (skipkey) continue;
   	   // we have a new key, in this and previous files: save it
	   ListAllKeys->Add(key);
         // read object from this source file
         next_source->cd( path );
         TObject *obj = key->ReadObj();
         if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
            // descendant of TH1 -> merge it
            //      cout << "Merging histogram " << obj->GetName() << endl;
//             TH1 *h1 = (TH1*)obj;
            // loop over all source files and add the content of the
            // correspondant histogram to the one pointed to by "h1"
         } else {
            // object is of no type that we know or can handle
            std::cout << "Unimplemented or unknown object type, name: "
            << obj->GetName() << " title: " << obj->GetTitle() << std::endl;
         }
         // now write the merged histogram (which is "in" obj) to the target file
         // note that this will just store obj in the current directory level,
         // which is not persistent until the complete directory itself is stored
         // by "target->Write()" below
         if ( obj ) {
            target->cd();
            //!!if the object is a tree, it is stored in globChain...
            if(obj->IsA()->InheritsFrom( TTree::Class() ))
               globChain->Merge(target->GetFile(),0,"keep");
            else
               obj->Write( key->GetName() );
         }
         oldkey=key;
      } // while ( ( TKey *key = (TKey*)nextkey() ) )
   } // while ( ( next_source = (TFile*)nextsource_iter() ) )
   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);
   std::cout << "Done." << std::endl;   
} 
