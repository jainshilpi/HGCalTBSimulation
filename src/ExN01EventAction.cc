#include "ExN01EventAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "ExN01CreateTree.hh"
#include "B1Run.hh"

#include <vector>

using namespace CLHEP;



EventAction::EventAction ()

: G4UserEventAction(),
  fEdep(0.)
 {} 



// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


EventAction::~EventAction ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void EventAction::BeginOfEventAction (const G4Event* evt)
{
  G4int evtNb = evt->GetEventID () ;
 std::cout<<"eventActio"<<std::endl;
  //CreateTree::Instance ()->Clear () ;
  // INSTANCE RUN/EVENT IN TREE
  //CreateTree::Instance ()->Event = evt->GetEventID () ;
 fEdep = 0.;

 //initialize tree variables to 0
 CreateTree::Instance ()->EnergyTotalAbs = 0 ;  

  ////put the energy in each layer separately
  for(int i=0; i<6; i++) ///make it better by making 4 a common variable
    {
      CreateTree::Instance ()->sensorE[i] = 0;
      //std::cout<<"At the bginning of event .... sensorE["<<i<<"] is "<<CreateTree::Instance ()->sensorE[i]<<std::endl;
    }

  
  for(int i=0; i<30; i++)
    {
      CreateTree::Instance ()->EnergyPerX0[i] = 0;
      CreateTree::Instance ()->depth_abs[i] = -9999;
    }

  ///

}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void EventAction::EndOfEventAction(const G4Event* evt)
{
float EnergyTotal;
  evt -> GetEventID();
  CreateTree::Instance ()->EnergyTotal = fEdep;
  B1Run* run  = static_cast<B1Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->AddEdep(fEdep);
//EnergyTotal = fEdep;

CreateTree::Instance ()->Fill ();
 std::cout<<EnergyTotal<<std::endl;
 //  EventAction::AddEdep(G4double edep);

 ////SJ printout

 for(int i=0; i<6; i++)
   std::cout<< CreateTree::Instance() -> sensorE[i] <<std::endl;
 
}

// ----------------




 





