#include "ExN01CreateTree.hh"
#include <algorithm>
#include <vector>

using namespace std ;

CreateTree* CreateTree::fInstance = NULL ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::CreateTree (TString name)
{
  if ( fInstance )
  {
    return ;
  }


  this->fInstance = this ;
  this->fname     = name ;
  
  this->ftree     = new TTree (name,name) ;
  
  this->ftree->SetDirectory(0);

  
  this->GetTree ()->Branch ("Event", &this->Event, "Event/I") ;
  
  this->GetTree ()->Branch ("EnergyTotalAbs", &this->EnergyTotalAbs, "EnergyTotalAbs/F") ;  



  //this->GetTree ()->Branch ("EnergyTotal", &this->EnergyTotal, "EnergyTotal/F") ;    

  this->GetTree ()->Branch ("sensorE", sensorE, "sensorE[6]/F") ;
  
  this->GetTree ()->Branch ("EnergyPerX0", EnergyPerX0, "EnergyPerX0[30]/F") ;
  this->GetTree ()->Branch ("depth_abs", depth_abs, "depth_abs[30]/F") ;  


  this->Clear () ;
  cout<<"tree"<<endl;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::~CreateTree ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int CreateTree::Fill () 
{ 
  return this->GetTree ()->Fill () ; 
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void CreateTree::Write (TFile * outfile)
{
  outfile->cd () ;
  ftree->Write () ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void CreateTree::Clear ()
{
  Event	= 0 ;

  EnergyTotalAbs = 0 ;  

  ////put the energy in each layer separately
  for(int i=0; i<6; i++) ///make it better by making 4 a common variable
    sensorE[i] = 0;


  
  for(int ii=0; ii<30; ii++)
    {
      EnergyPerX0[ii] = 0;
      depth_abs[ii] = -9999;
    }
  
 }
