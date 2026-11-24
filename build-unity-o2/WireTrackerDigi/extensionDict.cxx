// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME extensionDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiCollection.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiLocalCollection.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiV2Collection.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiAssociationCollection.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiV2AssociationCollection.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHitCollection.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHitSimTrackerHitLinkCollection.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/TrackCollection.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DatamodelDefinition.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigi.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiCollectionData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiLocal.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiLocalCollectionData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiLocalData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiLocalObj.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiObj.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiV2.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiV2CollectionData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiV2Data.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/DriftChamberDigiV2Obj.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiAssociation.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiAssociationCollectionData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiAssociationData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiAssociationObj.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiV2Association.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiV2AssociationCollectionData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiV2AssociationData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MCRecoDriftChamberDigiV2AssociationObj.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MutableDriftChamberDigi.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MutableDriftChamberDigiLocal.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MutableDriftChamberDigiV2.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MutableMCRecoDriftChamberDigiAssociation.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MutableMCRecoDriftChamberDigiV2Association.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MutableSenseWireHit.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MutableSenseWireHitSimTrackerHitLink.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/MutableTrack.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHit.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHitCollectionData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHitData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHitObj.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHitSimTrackerHitLink.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHitSimTrackerHitLinkCollectionData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHitSimTrackerHitLinkData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/SenseWireHitSimTrackerHitLinkObj.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/Track.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/TrackCollectionData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/TrackData.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/TrackObj.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/TrackerHit.h"
#include "/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/extension/TrackerHit_dev.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *extensioncLcLDriftChamberDigiData_Dictionary();
   static void extensioncLcLDriftChamberDigiData_TClassManip(TClass*);
   static void *new_extensioncLcLDriftChamberDigiData(void *p = nullptr);
   static void *newArray_extensioncLcLDriftChamberDigiData(Long_t size, void *p);
   static void delete_extensioncLcLDriftChamberDigiData(void *p);
   static void deleteArray_extensioncLcLDriftChamberDigiData(void *p);
   static void destruct_extensioncLcLDriftChamberDigiData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::DriftChamberDigiData*)
   {
      ::extension::DriftChamberDigiData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::DriftChamberDigiData));
      static ::ROOT::TGenericClassInfo 
         instance("extension::DriftChamberDigiData", 1, "extension/DriftChamberDigiData.h", 20,
                  typeid(::extension::DriftChamberDigiData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLDriftChamberDigiData_Dictionary, isa_proxy, 12,
                  sizeof(::extension::DriftChamberDigiData), alignof(::extension::DriftChamberDigiData) );
      instance.SetNew(&new_extensioncLcLDriftChamberDigiData);
      instance.SetNewArray(&newArray_extensioncLcLDriftChamberDigiData);
      instance.SetDelete(&delete_extensioncLcLDriftChamberDigiData);
      instance.SetDeleteArray(&deleteArray_extensioncLcLDriftChamberDigiData);
      instance.SetDestructor(&destruct_extensioncLcLDriftChamberDigiData);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("extension::DriftChamberDigiData","extension::v1::DriftChamberDigiData"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::DriftChamberDigiData*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::DriftChamberDigiData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLDriftChamberDigiData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiData*>(nullptr))->GetClass();
      extensioncLcLDriftChamberDigiData_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLDriftChamberDigiData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLDriftChamberDigi_Dictionary();
   static void extensioncLcLDriftChamberDigi_TClassManip(TClass*);
   static void *new_extensioncLcLDriftChamberDigi(void *p = nullptr);
   static void *newArray_extensioncLcLDriftChamberDigi(Long_t size, void *p);
   static void delete_extensioncLcLDriftChamberDigi(void *p);
   static void deleteArray_extensioncLcLDriftChamberDigi(void *p);
   static void destruct_extensioncLcLDriftChamberDigi(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::DriftChamberDigi*)
   {
      ::extension::DriftChamberDigi *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::DriftChamberDigi));
      static ::ROOT::TGenericClassInfo 
         instance("extension::DriftChamberDigi", 1, "extension/DriftChamberDigi.h", 45,
                  typeid(::extension::DriftChamberDigi), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLDriftChamberDigi_Dictionary, isa_proxy, 12,
                  sizeof(::extension::DriftChamberDigi), alignof(::extension::DriftChamberDigi) );
      instance.SetNew(&new_extensioncLcLDriftChamberDigi);
      instance.SetNewArray(&newArray_extensioncLcLDriftChamberDigi);
      instance.SetDelete(&delete_extensioncLcLDriftChamberDigi);
      instance.SetDeleteArray(&deleteArray_extensioncLcLDriftChamberDigi);
      instance.SetDestructor(&destruct_extensioncLcLDriftChamberDigi);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::DriftChamberDigi*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::DriftChamberDigi*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigi*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLDriftChamberDigi_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigi*>(nullptr))->GetClass();
      extensioncLcLDriftChamberDigi_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLDriftChamberDigi_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMutableDriftChamberDigi_Dictionary();
   static void extensioncLcLMutableDriftChamberDigi_TClassManip(TClass*);
   static void *new_extensioncLcLMutableDriftChamberDigi(void *p = nullptr);
   static void *newArray_extensioncLcLMutableDriftChamberDigi(Long_t size, void *p);
   static void delete_extensioncLcLMutableDriftChamberDigi(void *p);
   static void deleteArray_extensioncLcLMutableDriftChamberDigi(void *p);
   static void destruct_extensioncLcLMutableDriftChamberDigi(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MutableDriftChamberDigi*)
   {
      ::extension::MutableDriftChamberDigi *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MutableDriftChamberDigi));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MutableDriftChamberDigi", 1, "extension/MutableDriftChamberDigi.h", 37,
                  typeid(::extension::MutableDriftChamberDigi), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMutableDriftChamberDigi_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MutableDriftChamberDigi), alignof(::extension::MutableDriftChamberDigi) );
      instance.SetNew(&new_extensioncLcLMutableDriftChamberDigi);
      instance.SetNewArray(&newArray_extensioncLcLMutableDriftChamberDigi);
      instance.SetDelete(&delete_extensioncLcLMutableDriftChamberDigi);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMutableDriftChamberDigi);
      instance.SetDestructor(&destruct_extensioncLcLMutableDriftChamberDigi);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MutableDriftChamberDigi*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MutableDriftChamberDigi*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MutableDriftChamberDigi*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMutableDriftChamberDigi_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MutableDriftChamberDigi*>(nullptr))->GetClass();
      extensioncLcLMutableDriftChamberDigi_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMutableDriftChamberDigi_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLDriftChamberDigiCollection_Dictionary();
   static void extensioncLcLDriftChamberDigiCollection_TClassManip(TClass*);
   static void *new_extensioncLcLDriftChamberDigiCollection(void *p = nullptr);
   static void *newArray_extensioncLcLDriftChamberDigiCollection(Long_t size, void *p);
   static void delete_extensioncLcLDriftChamberDigiCollection(void *p);
   static void deleteArray_extensioncLcLDriftChamberDigiCollection(void *p);
   static void destruct_extensioncLcLDriftChamberDigiCollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::DriftChamberDigiCollection*)
   {
      ::extension::DriftChamberDigiCollection *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::DriftChamberDigiCollection));
      static ::ROOT::TGenericClassInfo 
         instance("extension::DriftChamberDigiCollection", 1, "", 148,
                  typeid(::extension::DriftChamberDigiCollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLDriftChamberDigiCollection_Dictionary, isa_proxy, 12,
                  sizeof(::extension::DriftChamberDigiCollection), alignof(::extension::DriftChamberDigiCollection) );
      instance.SetNew(&new_extensioncLcLDriftChamberDigiCollection);
      instance.SetNewArray(&newArray_extensioncLcLDriftChamberDigiCollection);
      instance.SetDelete(&delete_extensioncLcLDriftChamberDigiCollection);
      instance.SetDeleteArray(&deleteArray_extensioncLcLDriftChamberDigiCollection);
      instance.SetDestructor(&destruct_extensioncLcLDriftChamberDigiCollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::DriftChamberDigiCollection*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::DriftChamberDigiCollection*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiCollection*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLDriftChamberDigiCollection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiCollection*>(nullptr))->GetClass();
      extensioncLcLDriftChamberDigiCollection_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLDriftChamberDigiCollection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLDriftChamberDigiLocalData_Dictionary();
   static void extensioncLcLDriftChamberDigiLocalData_TClassManip(TClass*);
   static void *new_extensioncLcLDriftChamberDigiLocalData(void *p = nullptr);
   static void *newArray_extensioncLcLDriftChamberDigiLocalData(Long_t size, void *p);
   static void delete_extensioncLcLDriftChamberDigiLocalData(void *p);
   static void deleteArray_extensioncLcLDriftChamberDigiLocalData(void *p);
   static void destruct_extensioncLcLDriftChamberDigiLocalData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::DriftChamberDigiLocalData*)
   {
      ::extension::DriftChamberDigiLocalData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::DriftChamberDigiLocalData));
      static ::ROOT::TGenericClassInfo 
         instance("extension::DriftChamberDigiLocalData", 1, "extension/DriftChamberDigiLocalData.h", 19,
                  typeid(::extension::DriftChamberDigiLocalData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLDriftChamberDigiLocalData_Dictionary, isa_proxy, 12,
                  sizeof(::extension::DriftChamberDigiLocalData), alignof(::extension::DriftChamberDigiLocalData) );
      instance.SetNew(&new_extensioncLcLDriftChamberDigiLocalData);
      instance.SetNewArray(&newArray_extensioncLcLDriftChamberDigiLocalData);
      instance.SetDelete(&delete_extensioncLcLDriftChamberDigiLocalData);
      instance.SetDeleteArray(&deleteArray_extensioncLcLDriftChamberDigiLocalData);
      instance.SetDestructor(&destruct_extensioncLcLDriftChamberDigiLocalData);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("extension::DriftChamberDigiLocalData","extension::v1::DriftChamberDigiLocalData"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::DriftChamberDigiLocalData*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::DriftChamberDigiLocalData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiLocalData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLDriftChamberDigiLocalData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiLocalData*>(nullptr))->GetClass();
      extensioncLcLDriftChamberDigiLocalData_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLDriftChamberDigiLocalData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLDriftChamberDigiLocal_Dictionary();
   static void extensioncLcLDriftChamberDigiLocal_TClassManip(TClass*);
   static void *new_extensioncLcLDriftChamberDigiLocal(void *p = nullptr);
   static void *newArray_extensioncLcLDriftChamberDigiLocal(Long_t size, void *p);
   static void delete_extensioncLcLDriftChamberDigiLocal(void *p);
   static void deleteArray_extensioncLcLDriftChamberDigiLocal(void *p);
   static void destruct_extensioncLcLDriftChamberDigiLocal(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::DriftChamberDigiLocal*)
   {
      ::extension::DriftChamberDigiLocal *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::DriftChamberDigiLocal));
      static ::ROOT::TGenericClassInfo 
         instance("extension::DriftChamberDigiLocal", 1, "extension/DriftChamberDigiLocal.h", 43,
                  typeid(::extension::DriftChamberDigiLocal), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLDriftChamberDigiLocal_Dictionary, isa_proxy, 12,
                  sizeof(::extension::DriftChamberDigiLocal), alignof(::extension::DriftChamberDigiLocal) );
      instance.SetNew(&new_extensioncLcLDriftChamberDigiLocal);
      instance.SetNewArray(&newArray_extensioncLcLDriftChamberDigiLocal);
      instance.SetDelete(&delete_extensioncLcLDriftChamberDigiLocal);
      instance.SetDeleteArray(&deleteArray_extensioncLcLDriftChamberDigiLocal);
      instance.SetDestructor(&destruct_extensioncLcLDriftChamberDigiLocal);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::DriftChamberDigiLocal*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::DriftChamberDigiLocal*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiLocal*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLDriftChamberDigiLocal_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiLocal*>(nullptr))->GetClass();
      extensioncLcLDriftChamberDigiLocal_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLDriftChamberDigiLocal_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMutableDriftChamberDigiLocal_Dictionary();
   static void extensioncLcLMutableDriftChamberDigiLocal_TClassManip(TClass*);
   static void *new_extensioncLcLMutableDriftChamberDigiLocal(void *p = nullptr);
   static void *newArray_extensioncLcLMutableDriftChamberDigiLocal(Long_t size, void *p);
   static void delete_extensioncLcLMutableDriftChamberDigiLocal(void *p);
   static void deleteArray_extensioncLcLMutableDriftChamberDigiLocal(void *p);
   static void destruct_extensioncLcLMutableDriftChamberDigiLocal(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MutableDriftChamberDigiLocal*)
   {
      ::extension::MutableDriftChamberDigiLocal *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MutableDriftChamberDigiLocal));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MutableDriftChamberDigiLocal", 1, "extension/MutableDriftChamberDigiLocal.h", 35,
                  typeid(::extension::MutableDriftChamberDigiLocal), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMutableDriftChamberDigiLocal_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MutableDriftChamberDigiLocal), alignof(::extension::MutableDriftChamberDigiLocal) );
      instance.SetNew(&new_extensioncLcLMutableDriftChamberDigiLocal);
      instance.SetNewArray(&newArray_extensioncLcLMutableDriftChamberDigiLocal);
      instance.SetDelete(&delete_extensioncLcLMutableDriftChamberDigiLocal);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMutableDriftChamberDigiLocal);
      instance.SetDestructor(&destruct_extensioncLcLMutableDriftChamberDigiLocal);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MutableDriftChamberDigiLocal*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MutableDriftChamberDigiLocal*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MutableDriftChamberDigiLocal*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMutableDriftChamberDigiLocal_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MutableDriftChamberDigiLocal*>(nullptr))->GetClass();
      extensioncLcLMutableDriftChamberDigiLocal_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMutableDriftChamberDigiLocal_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLDriftChamberDigiLocalCollection_Dictionary();
   static void extensioncLcLDriftChamberDigiLocalCollection_TClassManip(TClass*);
   static void *new_extensioncLcLDriftChamberDigiLocalCollection(void *p = nullptr);
   static void *newArray_extensioncLcLDriftChamberDigiLocalCollection(Long_t size, void *p);
   static void delete_extensioncLcLDriftChamberDigiLocalCollection(void *p);
   static void deleteArray_extensioncLcLDriftChamberDigiLocalCollection(void *p);
   static void destruct_extensioncLcLDriftChamberDigiLocalCollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::DriftChamberDigiLocalCollection*)
   {
      ::extension::DriftChamberDigiLocalCollection *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::DriftChamberDigiLocalCollection));
      static ::ROOT::TGenericClassInfo 
         instance("extension::DriftChamberDigiLocalCollection", 1, "", 552,
                  typeid(::extension::DriftChamberDigiLocalCollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLDriftChamberDigiLocalCollection_Dictionary, isa_proxy, 12,
                  sizeof(::extension::DriftChamberDigiLocalCollection), alignof(::extension::DriftChamberDigiLocalCollection) );
      instance.SetNew(&new_extensioncLcLDriftChamberDigiLocalCollection);
      instance.SetNewArray(&newArray_extensioncLcLDriftChamberDigiLocalCollection);
      instance.SetDelete(&delete_extensioncLcLDriftChamberDigiLocalCollection);
      instance.SetDeleteArray(&deleteArray_extensioncLcLDriftChamberDigiLocalCollection);
      instance.SetDestructor(&destruct_extensioncLcLDriftChamberDigiLocalCollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::DriftChamberDigiLocalCollection*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::DriftChamberDigiLocalCollection*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiLocalCollection*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLDriftChamberDigiLocalCollection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiLocalCollection*>(nullptr))->GetClass();
      extensioncLcLDriftChamberDigiLocalCollection_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLDriftChamberDigiLocalCollection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLDriftChamberDigiV2Data_Dictionary();
   static void extensioncLcLDriftChamberDigiV2Data_TClassManip(TClass*);
   static void *new_extensioncLcLDriftChamberDigiV2Data(void *p = nullptr);
   static void *newArray_extensioncLcLDriftChamberDigiV2Data(Long_t size, void *p);
   static void delete_extensioncLcLDriftChamberDigiV2Data(void *p);
   static void deleteArray_extensioncLcLDriftChamberDigiV2Data(void *p);
   static void destruct_extensioncLcLDriftChamberDigiV2Data(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::DriftChamberDigiV2Data*)
   {
      ::extension::DriftChamberDigiV2Data *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::DriftChamberDigiV2Data));
      static ::ROOT::TGenericClassInfo 
         instance("extension::DriftChamberDigiV2Data", 1, "extension/DriftChamberDigiV2Data.h", 20,
                  typeid(::extension::DriftChamberDigiV2Data), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLDriftChamberDigiV2Data_Dictionary, isa_proxy, 12,
                  sizeof(::extension::DriftChamberDigiV2Data), alignof(::extension::DriftChamberDigiV2Data) );
      instance.SetNew(&new_extensioncLcLDriftChamberDigiV2Data);
      instance.SetNewArray(&newArray_extensioncLcLDriftChamberDigiV2Data);
      instance.SetDelete(&delete_extensioncLcLDriftChamberDigiV2Data);
      instance.SetDeleteArray(&deleteArray_extensioncLcLDriftChamberDigiV2Data);
      instance.SetDestructor(&destruct_extensioncLcLDriftChamberDigiV2Data);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("extension::DriftChamberDigiV2Data","extension::v1::DriftChamberDigiV2Data"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::DriftChamberDigiV2Data*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::DriftChamberDigiV2Data*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiV2Data*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLDriftChamberDigiV2Data_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiV2Data*>(nullptr))->GetClass();
      extensioncLcLDriftChamberDigiV2Data_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLDriftChamberDigiV2Data_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLDriftChamberDigiV2_Dictionary();
   static void extensioncLcLDriftChamberDigiV2_TClassManip(TClass*);
   static void *new_extensioncLcLDriftChamberDigiV2(void *p = nullptr);
   static void *newArray_extensioncLcLDriftChamberDigiV2(Long_t size, void *p);
   static void delete_extensioncLcLDriftChamberDigiV2(void *p);
   static void deleteArray_extensioncLcLDriftChamberDigiV2(void *p);
   static void destruct_extensioncLcLDriftChamberDigiV2(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::DriftChamberDigiV2*)
   {
      ::extension::DriftChamberDigiV2 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::DriftChamberDigiV2));
      static ::ROOT::TGenericClassInfo 
         instance("extension::DriftChamberDigiV2", 1, "extension/DriftChamberDigiV2.h", 46,
                  typeid(::extension::DriftChamberDigiV2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLDriftChamberDigiV2_Dictionary, isa_proxy, 12,
                  sizeof(::extension::DriftChamberDigiV2), alignof(::extension::DriftChamberDigiV2) );
      instance.SetNew(&new_extensioncLcLDriftChamberDigiV2);
      instance.SetNewArray(&newArray_extensioncLcLDriftChamberDigiV2);
      instance.SetDelete(&delete_extensioncLcLDriftChamberDigiV2);
      instance.SetDeleteArray(&deleteArray_extensioncLcLDriftChamberDigiV2);
      instance.SetDestructor(&destruct_extensioncLcLDriftChamberDigiV2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::DriftChamberDigiV2*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::DriftChamberDigiV2*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiV2*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLDriftChamberDigiV2_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiV2*>(nullptr))->GetClass();
      extensioncLcLDriftChamberDigiV2_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLDriftChamberDigiV2_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMutableDriftChamberDigiV2_Dictionary();
   static void extensioncLcLMutableDriftChamberDigiV2_TClassManip(TClass*);
   static void *new_extensioncLcLMutableDriftChamberDigiV2(void *p = nullptr);
   static void *newArray_extensioncLcLMutableDriftChamberDigiV2(Long_t size, void *p);
   static void delete_extensioncLcLMutableDriftChamberDigiV2(void *p);
   static void deleteArray_extensioncLcLMutableDriftChamberDigiV2(void *p);
   static void destruct_extensioncLcLMutableDriftChamberDigiV2(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MutableDriftChamberDigiV2*)
   {
      ::extension::MutableDriftChamberDigiV2 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MutableDriftChamberDigiV2));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MutableDriftChamberDigiV2", 1, "extension/MutableDriftChamberDigiV2.h", 38,
                  typeid(::extension::MutableDriftChamberDigiV2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMutableDriftChamberDigiV2_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MutableDriftChamberDigiV2), alignof(::extension::MutableDriftChamberDigiV2) );
      instance.SetNew(&new_extensioncLcLMutableDriftChamberDigiV2);
      instance.SetNewArray(&newArray_extensioncLcLMutableDriftChamberDigiV2);
      instance.SetDelete(&delete_extensioncLcLMutableDriftChamberDigiV2);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMutableDriftChamberDigiV2);
      instance.SetDestructor(&destruct_extensioncLcLMutableDriftChamberDigiV2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MutableDriftChamberDigiV2*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MutableDriftChamberDigiV2*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MutableDriftChamberDigiV2*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMutableDriftChamberDigiV2_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MutableDriftChamberDigiV2*>(nullptr))->GetClass();
      extensioncLcLMutableDriftChamberDigiV2_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMutableDriftChamberDigiV2_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLDriftChamberDigiV2Collection_Dictionary();
   static void extensioncLcLDriftChamberDigiV2Collection_TClassManip(TClass*);
   static void *new_extensioncLcLDriftChamberDigiV2Collection(void *p = nullptr);
   static void *newArray_extensioncLcLDriftChamberDigiV2Collection(Long_t size, void *p);
   static void delete_extensioncLcLDriftChamberDigiV2Collection(void *p);
   static void deleteArray_extensioncLcLDriftChamberDigiV2Collection(void *p);
   static void destruct_extensioncLcLDriftChamberDigiV2Collection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::DriftChamberDigiV2Collection*)
   {
      ::extension::DriftChamberDigiV2Collection *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::DriftChamberDigiV2Collection));
      static ::ROOT::TGenericClassInfo 
         instance("extension::DriftChamberDigiV2Collection", 1, "", 953,
                  typeid(::extension::DriftChamberDigiV2Collection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLDriftChamberDigiV2Collection_Dictionary, isa_proxy, 12,
                  sizeof(::extension::DriftChamberDigiV2Collection), alignof(::extension::DriftChamberDigiV2Collection) );
      instance.SetNew(&new_extensioncLcLDriftChamberDigiV2Collection);
      instance.SetNewArray(&newArray_extensioncLcLDriftChamberDigiV2Collection);
      instance.SetDelete(&delete_extensioncLcLDriftChamberDigiV2Collection);
      instance.SetDeleteArray(&deleteArray_extensioncLcLDriftChamberDigiV2Collection);
      instance.SetDestructor(&destruct_extensioncLcLDriftChamberDigiV2Collection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::DriftChamberDigiV2Collection*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::DriftChamberDigiV2Collection*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiV2Collection*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLDriftChamberDigiV2Collection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::DriftChamberDigiV2Collection*>(nullptr))->GetClass();
      extensioncLcLDriftChamberDigiV2Collection_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLDriftChamberDigiV2Collection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMCRecoDriftChamberDigiAssociationData_Dictionary();
   static void extensioncLcLMCRecoDriftChamberDigiAssociationData_TClassManip(TClass*);
   static void *new_extensioncLcLMCRecoDriftChamberDigiAssociationData(void *p = nullptr);
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiAssociationData(Long_t size, void *p);
   static void delete_extensioncLcLMCRecoDriftChamberDigiAssociationData(void *p);
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiAssociationData(void *p);
   static void destruct_extensioncLcLMCRecoDriftChamberDigiAssociationData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MCRecoDriftChamberDigiAssociationData*)
   {
      ::extension::MCRecoDriftChamberDigiAssociationData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MCRecoDriftChamberDigiAssociationData));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MCRecoDriftChamberDigiAssociationData", 1, "extension/MCRecoDriftChamberDigiAssociationData.h", 18,
                  typeid(::extension::MCRecoDriftChamberDigiAssociationData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMCRecoDriftChamberDigiAssociationData_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MCRecoDriftChamberDigiAssociationData), alignof(::extension::MCRecoDriftChamberDigiAssociationData) );
      instance.SetNew(&new_extensioncLcLMCRecoDriftChamberDigiAssociationData);
      instance.SetNewArray(&newArray_extensioncLcLMCRecoDriftChamberDigiAssociationData);
      instance.SetDelete(&delete_extensioncLcLMCRecoDriftChamberDigiAssociationData);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMCRecoDriftChamberDigiAssociationData);
      instance.SetDestructor(&destruct_extensioncLcLMCRecoDriftChamberDigiAssociationData);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("extension::MCRecoDriftChamberDigiAssociationData","extension::v1::MCRecoDriftChamberDigiAssociationData"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MCRecoDriftChamberDigiAssociationData*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MCRecoDriftChamberDigiAssociationData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiAssociationData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMCRecoDriftChamberDigiAssociationData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiAssociationData*>(nullptr))->GetClass();
      extensioncLcLMCRecoDriftChamberDigiAssociationData_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMCRecoDriftChamberDigiAssociationData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMCRecoDriftChamberDigiAssociation_Dictionary();
   static void extensioncLcLMCRecoDriftChamberDigiAssociation_TClassManip(TClass*);
   static void *new_extensioncLcLMCRecoDriftChamberDigiAssociation(void *p = nullptr);
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiAssociation(Long_t size, void *p);
   static void delete_extensioncLcLMCRecoDriftChamberDigiAssociation(void *p);
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiAssociation(void *p);
   static void destruct_extensioncLcLMCRecoDriftChamberDigiAssociation(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MCRecoDriftChamberDigiAssociation*)
   {
      ::extension::MCRecoDriftChamberDigiAssociation *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MCRecoDriftChamberDigiAssociation));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MCRecoDriftChamberDigiAssociation", 1, "extension/MCRecoDriftChamberDigiAssociation.h", 48,
                  typeid(::extension::MCRecoDriftChamberDigiAssociation), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMCRecoDriftChamberDigiAssociation_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MCRecoDriftChamberDigiAssociation), alignof(::extension::MCRecoDriftChamberDigiAssociation) );
      instance.SetNew(&new_extensioncLcLMCRecoDriftChamberDigiAssociation);
      instance.SetNewArray(&newArray_extensioncLcLMCRecoDriftChamberDigiAssociation);
      instance.SetDelete(&delete_extensioncLcLMCRecoDriftChamberDigiAssociation);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMCRecoDriftChamberDigiAssociation);
      instance.SetDestructor(&destruct_extensioncLcLMCRecoDriftChamberDigiAssociation);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MCRecoDriftChamberDigiAssociation*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MCRecoDriftChamberDigiAssociation*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiAssociation*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMCRecoDriftChamberDigiAssociation_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiAssociation*>(nullptr))->GetClass();
      extensioncLcLMCRecoDriftChamberDigiAssociation_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMCRecoDriftChamberDigiAssociation_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMutableMCRecoDriftChamberDigiAssociation_Dictionary();
   static void extensioncLcLMutableMCRecoDriftChamberDigiAssociation_TClassManip(TClass*);
   static void *new_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(void *p = nullptr);
   static void *newArray_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(Long_t size, void *p);
   static void delete_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(void *p);
   static void deleteArray_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(void *p);
   static void destruct_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MutableMCRecoDriftChamberDigiAssociation*)
   {
      ::extension::MutableMCRecoDriftChamberDigiAssociation *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MutableMCRecoDriftChamberDigiAssociation));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MutableMCRecoDriftChamberDigiAssociation", 1, "extension/MutableMCRecoDriftChamberDigiAssociation.h", 40,
                  typeid(::extension::MutableMCRecoDriftChamberDigiAssociation), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMutableMCRecoDriftChamberDigiAssociation_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MutableMCRecoDriftChamberDigiAssociation), alignof(::extension::MutableMCRecoDriftChamberDigiAssociation) );
      instance.SetNew(&new_extensioncLcLMutableMCRecoDriftChamberDigiAssociation);
      instance.SetNewArray(&newArray_extensioncLcLMutableMCRecoDriftChamberDigiAssociation);
      instance.SetDelete(&delete_extensioncLcLMutableMCRecoDriftChamberDigiAssociation);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMutableMCRecoDriftChamberDigiAssociation);
      instance.SetDestructor(&destruct_extensioncLcLMutableMCRecoDriftChamberDigiAssociation);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MutableMCRecoDriftChamberDigiAssociation*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MutableMCRecoDriftChamberDigiAssociation*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MutableMCRecoDriftChamberDigiAssociation*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMutableMCRecoDriftChamberDigiAssociation_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MutableMCRecoDriftChamberDigiAssociation*>(nullptr))->GetClass();
      extensioncLcLMutableMCRecoDriftChamberDigiAssociation_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMutableMCRecoDriftChamberDigiAssociation_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMCRecoDriftChamberDigiAssociationCollection_Dictionary();
   static void extensioncLcLMCRecoDriftChamberDigiAssociationCollection_TClassManip(TClass*);
   static void *new_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(void *p = nullptr);
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(Long_t size, void *p);
   static void delete_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(void *p);
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(void *p);
   static void destruct_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MCRecoDriftChamberDigiAssociationCollection*)
   {
      ::extension::MCRecoDriftChamberDigiAssociationCollection *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MCRecoDriftChamberDigiAssociationCollection));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MCRecoDriftChamberDigiAssociationCollection", 1, "", 1360,
                  typeid(::extension::MCRecoDriftChamberDigiAssociationCollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMCRecoDriftChamberDigiAssociationCollection_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MCRecoDriftChamberDigiAssociationCollection), alignof(::extension::MCRecoDriftChamberDigiAssociationCollection) );
      instance.SetNew(&new_extensioncLcLMCRecoDriftChamberDigiAssociationCollection);
      instance.SetNewArray(&newArray_extensioncLcLMCRecoDriftChamberDigiAssociationCollection);
      instance.SetDelete(&delete_extensioncLcLMCRecoDriftChamberDigiAssociationCollection);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMCRecoDriftChamberDigiAssociationCollection);
      instance.SetDestructor(&destruct_extensioncLcLMCRecoDriftChamberDigiAssociationCollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MCRecoDriftChamberDigiAssociationCollection*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MCRecoDriftChamberDigiAssociationCollection*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiAssociationCollection*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMCRecoDriftChamberDigiAssociationCollection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiAssociationCollection*>(nullptr))->GetClass();
      extensioncLcLMCRecoDriftChamberDigiAssociationCollection_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMCRecoDriftChamberDigiAssociationCollection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMCRecoDriftChamberDigiV2AssociationData_Dictionary();
   static void extensioncLcLMCRecoDriftChamberDigiV2AssociationData_TClassManip(TClass*);
   static void *new_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(void *p = nullptr);
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(Long_t size, void *p);
   static void delete_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(void *p);
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(void *p);
   static void destruct_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MCRecoDriftChamberDigiV2AssociationData*)
   {
      ::extension::MCRecoDriftChamberDigiV2AssociationData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MCRecoDriftChamberDigiV2AssociationData));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MCRecoDriftChamberDigiV2AssociationData", 1, "extension/MCRecoDriftChamberDigiV2AssociationData.h", 18,
                  typeid(::extension::MCRecoDriftChamberDigiV2AssociationData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMCRecoDriftChamberDigiV2AssociationData_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MCRecoDriftChamberDigiV2AssociationData), alignof(::extension::MCRecoDriftChamberDigiV2AssociationData) );
      instance.SetNew(&new_extensioncLcLMCRecoDriftChamberDigiV2AssociationData);
      instance.SetNewArray(&newArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationData);
      instance.SetDelete(&delete_extensioncLcLMCRecoDriftChamberDigiV2AssociationData);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationData);
      instance.SetDestructor(&destruct_extensioncLcLMCRecoDriftChamberDigiV2AssociationData);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("extension::MCRecoDriftChamberDigiV2AssociationData","extension::v1::MCRecoDriftChamberDigiV2AssociationData"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MCRecoDriftChamberDigiV2AssociationData*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MCRecoDriftChamberDigiV2AssociationData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiV2AssociationData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMCRecoDriftChamberDigiV2AssociationData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiV2AssociationData*>(nullptr))->GetClass();
      extensioncLcLMCRecoDriftChamberDigiV2AssociationData_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMCRecoDriftChamberDigiV2AssociationData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMCRecoDriftChamberDigiV2Association_Dictionary();
   static void extensioncLcLMCRecoDriftChamberDigiV2Association_TClassManip(TClass*);
   static void *new_extensioncLcLMCRecoDriftChamberDigiV2Association(void *p = nullptr);
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiV2Association(Long_t size, void *p);
   static void delete_extensioncLcLMCRecoDriftChamberDigiV2Association(void *p);
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiV2Association(void *p);
   static void destruct_extensioncLcLMCRecoDriftChamberDigiV2Association(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MCRecoDriftChamberDigiV2Association*)
   {
      ::extension::MCRecoDriftChamberDigiV2Association *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MCRecoDriftChamberDigiV2Association));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MCRecoDriftChamberDigiV2Association", 1, "extension/MCRecoDriftChamberDigiV2Association.h", 48,
                  typeid(::extension::MCRecoDriftChamberDigiV2Association), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMCRecoDriftChamberDigiV2Association_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MCRecoDriftChamberDigiV2Association), alignof(::extension::MCRecoDriftChamberDigiV2Association) );
      instance.SetNew(&new_extensioncLcLMCRecoDriftChamberDigiV2Association);
      instance.SetNewArray(&newArray_extensioncLcLMCRecoDriftChamberDigiV2Association);
      instance.SetDelete(&delete_extensioncLcLMCRecoDriftChamberDigiV2Association);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMCRecoDriftChamberDigiV2Association);
      instance.SetDestructor(&destruct_extensioncLcLMCRecoDriftChamberDigiV2Association);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MCRecoDriftChamberDigiV2Association*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MCRecoDriftChamberDigiV2Association*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiV2Association*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMCRecoDriftChamberDigiV2Association_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiV2Association*>(nullptr))->GetClass();
      extensioncLcLMCRecoDriftChamberDigiV2Association_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMCRecoDriftChamberDigiV2Association_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMutableMCRecoDriftChamberDigiV2Association_Dictionary();
   static void extensioncLcLMutableMCRecoDriftChamberDigiV2Association_TClassManip(TClass*);
   static void *new_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(void *p = nullptr);
   static void *newArray_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(Long_t size, void *p);
   static void delete_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(void *p);
   static void deleteArray_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(void *p);
   static void destruct_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MutableMCRecoDriftChamberDigiV2Association*)
   {
      ::extension::MutableMCRecoDriftChamberDigiV2Association *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MutableMCRecoDriftChamberDigiV2Association));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MutableMCRecoDriftChamberDigiV2Association", 1, "extension/MutableMCRecoDriftChamberDigiV2Association.h", 40,
                  typeid(::extension::MutableMCRecoDriftChamberDigiV2Association), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMutableMCRecoDriftChamberDigiV2Association_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MutableMCRecoDriftChamberDigiV2Association), alignof(::extension::MutableMCRecoDriftChamberDigiV2Association) );
      instance.SetNew(&new_extensioncLcLMutableMCRecoDriftChamberDigiV2Association);
      instance.SetNewArray(&newArray_extensioncLcLMutableMCRecoDriftChamberDigiV2Association);
      instance.SetDelete(&delete_extensioncLcLMutableMCRecoDriftChamberDigiV2Association);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMutableMCRecoDriftChamberDigiV2Association);
      instance.SetDestructor(&destruct_extensioncLcLMutableMCRecoDriftChamberDigiV2Association);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MutableMCRecoDriftChamberDigiV2Association*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MutableMCRecoDriftChamberDigiV2Association*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MutableMCRecoDriftChamberDigiV2Association*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMutableMCRecoDriftChamberDigiV2Association_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MutableMCRecoDriftChamberDigiV2Association*>(nullptr))->GetClass();
      extensioncLcLMutableMCRecoDriftChamberDigiV2Association_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMutableMCRecoDriftChamberDigiV2Association_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection_Dictionary();
   static void extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection_TClassManip(TClass*);
   static void *new_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(void *p = nullptr);
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(Long_t size, void *p);
   static void delete_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(void *p);
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(void *p);
   static void destruct_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MCRecoDriftChamberDigiV2AssociationCollection*)
   {
      ::extension::MCRecoDriftChamberDigiV2AssociationCollection *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MCRecoDriftChamberDigiV2AssociationCollection));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MCRecoDriftChamberDigiV2AssociationCollection", 1, "", 1755,
                  typeid(::extension::MCRecoDriftChamberDigiV2AssociationCollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MCRecoDriftChamberDigiV2AssociationCollection), alignof(::extension::MCRecoDriftChamberDigiV2AssociationCollection) );
      instance.SetNew(&new_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection);
      instance.SetNewArray(&newArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection);
      instance.SetDelete(&delete_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection);
      instance.SetDestructor(&destruct_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MCRecoDriftChamberDigiV2AssociationCollection*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MCRecoDriftChamberDigiV2AssociationCollection*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiV2AssociationCollection*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MCRecoDriftChamberDigiV2AssociationCollection*>(nullptr))->GetClass();
      extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLSenseWireHitData_Dictionary();
   static void extensioncLcLSenseWireHitData_TClassManip(TClass*);
   static void *new_extensioncLcLSenseWireHitData(void *p = nullptr);
   static void *newArray_extensioncLcLSenseWireHitData(Long_t size, void *p);
   static void delete_extensioncLcLSenseWireHitData(void *p);
   static void deleteArray_extensioncLcLSenseWireHitData(void *p);
   static void destruct_extensioncLcLSenseWireHitData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::SenseWireHitData*)
   {
      ::extension::SenseWireHitData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::SenseWireHitData));
      static ::ROOT::TGenericClassInfo 
         instance("extension::SenseWireHitData", 1, "extension/SenseWireHitData.h", 20,
                  typeid(::extension::SenseWireHitData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLSenseWireHitData_Dictionary, isa_proxy, 12,
                  sizeof(::extension::SenseWireHitData), alignof(::extension::SenseWireHitData) );
      instance.SetNew(&new_extensioncLcLSenseWireHitData);
      instance.SetNewArray(&newArray_extensioncLcLSenseWireHitData);
      instance.SetDelete(&delete_extensioncLcLSenseWireHitData);
      instance.SetDeleteArray(&deleteArray_extensioncLcLSenseWireHitData);
      instance.SetDestructor(&destruct_extensioncLcLSenseWireHitData);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("extension::SenseWireHitData","extension::v1::SenseWireHitData"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::SenseWireHitData*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::SenseWireHitData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLSenseWireHitData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitData*>(nullptr))->GetClass();
      extensioncLcLSenseWireHitData_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLSenseWireHitData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLSenseWireHit_Dictionary();
   static void extensioncLcLSenseWireHit_TClassManip(TClass*);
   static void *new_extensioncLcLSenseWireHit(void *p = nullptr);
   static void *newArray_extensioncLcLSenseWireHit(Long_t size, void *p);
   static void delete_extensioncLcLSenseWireHit(void *p);
   static void deleteArray_extensioncLcLSenseWireHit(void *p);
   static void destruct_extensioncLcLSenseWireHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::SenseWireHit*)
   {
      ::extension::SenseWireHit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::SenseWireHit));
      static ::ROOT::TGenericClassInfo 
         instance("extension::SenseWireHit", 1, "extension/SenseWireHit.h", 48,
                  typeid(::extension::SenseWireHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLSenseWireHit_Dictionary, isa_proxy, 12,
                  sizeof(::extension::SenseWireHit), alignof(::extension::SenseWireHit) );
      instance.SetNew(&new_extensioncLcLSenseWireHit);
      instance.SetNewArray(&newArray_extensioncLcLSenseWireHit);
      instance.SetDelete(&delete_extensioncLcLSenseWireHit);
      instance.SetDeleteArray(&deleteArray_extensioncLcLSenseWireHit);
      instance.SetDestructor(&destruct_extensioncLcLSenseWireHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::SenseWireHit*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::SenseWireHit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLSenseWireHit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHit*>(nullptr))->GetClass();
      extensioncLcLSenseWireHit_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLSenseWireHit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMutableSenseWireHit_Dictionary();
   static void extensioncLcLMutableSenseWireHit_TClassManip(TClass*);
   static void *new_extensioncLcLMutableSenseWireHit(void *p = nullptr);
   static void *newArray_extensioncLcLMutableSenseWireHit(Long_t size, void *p);
   static void delete_extensioncLcLMutableSenseWireHit(void *p);
   static void deleteArray_extensioncLcLMutableSenseWireHit(void *p);
   static void destruct_extensioncLcLMutableSenseWireHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MutableSenseWireHit*)
   {
      ::extension::MutableSenseWireHit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MutableSenseWireHit));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MutableSenseWireHit", 1, "extension/MutableSenseWireHit.h", 40,
                  typeid(::extension::MutableSenseWireHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMutableSenseWireHit_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MutableSenseWireHit), alignof(::extension::MutableSenseWireHit) );
      instance.SetNew(&new_extensioncLcLMutableSenseWireHit);
      instance.SetNewArray(&newArray_extensioncLcLMutableSenseWireHit);
      instance.SetDelete(&delete_extensioncLcLMutableSenseWireHit);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMutableSenseWireHit);
      instance.SetDestructor(&destruct_extensioncLcLMutableSenseWireHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MutableSenseWireHit*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MutableSenseWireHit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MutableSenseWireHit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMutableSenseWireHit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MutableSenseWireHit*>(nullptr))->GetClass();
      extensioncLcLMutableSenseWireHit_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMutableSenseWireHit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLSenseWireHitCollection_Dictionary();
   static void extensioncLcLSenseWireHitCollection_TClassManip(TClass*);
   static void *new_extensioncLcLSenseWireHitCollection(void *p = nullptr);
   static void *newArray_extensioncLcLSenseWireHitCollection(Long_t size, void *p);
   static void delete_extensioncLcLSenseWireHitCollection(void *p);
   static void deleteArray_extensioncLcLSenseWireHitCollection(void *p);
   static void destruct_extensioncLcLSenseWireHitCollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::SenseWireHitCollection*)
   {
      ::extension::SenseWireHitCollection *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::SenseWireHitCollection));
      static ::ROOT::TGenericClassInfo 
         instance("extension::SenseWireHitCollection", 1, "", 2150,
                  typeid(::extension::SenseWireHitCollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLSenseWireHitCollection_Dictionary, isa_proxy, 12,
                  sizeof(::extension::SenseWireHitCollection), alignof(::extension::SenseWireHitCollection) );
      instance.SetNew(&new_extensioncLcLSenseWireHitCollection);
      instance.SetNewArray(&newArray_extensioncLcLSenseWireHitCollection);
      instance.SetDelete(&delete_extensioncLcLSenseWireHitCollection);
      instance.SetDeleteArray(&deleteArray_extensioncLcLSenseWireHitCollection);
      instance.SetDestructor(&destruct_extensioncLcLSenseWireHitCollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::SenseWireHitCollection*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::SenseWireHitCollection*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitCollection*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLSenseWireHitCollection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitCollection*>(nullptr))->GetClass();
      extensioncLcLSenseWireHitCollection_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLSenseWireHitCollection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLSenseWireHitSimTrackerHitLinkData_Dictionary();
   static void extensioncLcLSenseWireHitSimTrackerHitLinkData_TClassManip(TClass*);
   static void *new_extensioncLcLSenseWireHitSimTrackerHitLinkData(void *p = nullptr);
   static void *newArray_extensioncLcLSenseWireHitSimTrackerHitLinkData(Long_t size, void *p);
   static void delete_extensioncLcLSenseWireHitSimTrackerHitLinkData(void *p);
   static void deleteArray_extensioncLcLSenseWireHitSimTrackerHitLinkData(void *p);
   static void destruct_extensioncLcLSenseWireHitSimTrackerHitLinkData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::SenseWireHitSimTrackerHitLinkData*)
   {
      ::extension::SenseWireHitSimTrackerHitLinkData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::SenseWireHitSimTrackerHitLinkData));
      static ::ROOT::TGenericClassInfo 
         instance("extension::SenseWireHitSimTrackerHitLinkData", 1, "extension/SenseWireHitSimTrackerHitLinkData.h", 18,
                  typeid(::extension::SenseWireHitSimTrackerHitLinkData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLSenseWireHitSimTrackerHitLinkData_Dictionary, isa_proxy, 12,
                  sizeof(::extension::SenseWireHitSimTrackerHitLinkData), alignof(::extension::SenseWireHitSimTrackerHitLinkData) );
      instance.SetNew(&new_extensioncLcLSenseWireHitSimTrackerHitLinkData);
      instance.SetNewArray(&newArray_extensioncLcLSenseWireHitSimTrackerHitLinkData);
      instance.SetDelete(&delete_extensioncLcLSenseWireHitSimTrackerHitLinkData);
      instance.SetDeleteArray(&deleteArray_extensioncLcLSenseWireHitSimTrackerHitLinkData);
      instance.SetDestructor(&destruct_extensioncLcLSenseWireHitSimTrackerHitLinkData);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("extension::SenseWireHitSimTrackerHitLinkData","extension::v1::SenseWireHitSimTrackerHitLinkData"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::SenseWireHitSimTrackerHitLinkData*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::SenseWireHitSimTrackerHitLinkData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitSimTrackerHitLinkData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLSenseWireHitSimTrackerHitLinkData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitSimTrackerHitLinkData*>(nullptr))->GetClass();
      extensioncLcLSenseWireHitSimTrackerHitLinkData_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLSenseWireHitSimTrackerHitLinkData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLSenseWireHitSimTrackerHitLink_Dictionary();
   static void extensioncLcLSenseWireHitSimTrackerHitLink_TClassManip(TClass*);
   static void *new_extensioncLcLSenseWireHitSimTrackerHitLink(void *p = nullptr);
   static void *newArray_extensioncLcLSenseWireHitSimTrackerHitLink(Long_t size, void *p);
   static void delete_extensioncLcLSenseWireHitSimTrackerHitLink(void *p);
   static void deleteArray_extensioncLcLSenseWireHitSimTrackerHitLink(void *p);
   static void destruct_extensioncLcLSenseWireHitSimTrackerHitLink(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::SenseWireHitSimTrackerHitLink*)
   {
      ::extension::SenseWireHitSimTrackerHitLink *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::SenseWireHitSimTrackerHitLink));
      static ::ROOT::TGenericClassInfo 
         instance("extension::SenseWireHitSimTrackerHitLink", 1, "extension/SenseWireHitSimTrackerHitLink.h", 48,
                  typeid(::extension::SenseWireHitSimTrackerHitLink), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLSenseWireHitSimTrackerHitLink_Dictionary, isa_proxy, 12,
                  sizeof(::extension::SenseWireHitSimTrackerHitLink), alignof(::extension::SenseWireHitSimTrackerHitLink) );
      instance.SetNew(&new_extensioncLcLSenseWireHitSimTrackerHitLink);
      instance.SetNewArray(&newArray_extensioncLcLSenseWireHitSimTrackerHitLink);
      instance.SetDelete(&delete_extensioncLcLSenseWireHitSimTrackerHitLink);
      instance.SetDeleteArray(&deleteArray_extensioncLcLSenseWireHitSimTrackerHitLink);
      instance.SetDestructor(&destruct_extensioncLcLSenseWireHitSimTrackerHitLink);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::SenseWireHitSimTrackerHitLink*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::SenseWireHitSimTrackerHitLink*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitSimTrackerHitLink*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLSenseWireHitSimTrackerHitLink_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitSimTrackerHitLink*>(nullptr))->GetClass();
      extensioncLcLSenseWireHitSimTrackerHitLink_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLSenseWireHitSimTrackerHitLink_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMutableSenseWireHitSimTrackerHitLink_Dictionary();
   static void extensioncLcLMutableSenseWireHitSimTrackerHitLink_TClassManip(TClass*);
   static void *new_extensioncLcLMutableSenseWireHitSimTrackerHitLink(void *p = nullptr);
   static void *newArray_extensioncLcLMutableSenseWireHitSimTrackerHitLink(Long_t size, void *p);
   static void delete_extensioncLcLMutableSenseWireHitSimTrackerHitLink(void *p);
   static void deleteArray_extensioncLcLMutableSenseWireHitSimTrackerHitLink(void *p);
   static void destruct_extensioncLcLMutableSenseWireHitSimTrackerHitLink(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MutableSenseWireHitSimTrackerHitLink*)
   {
      ::extension::MutableSenseWireHitSimTrackerHitLink *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MutableSenseWireHitSimTrackerHitLink));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MutableSenseWireHitSimTrackerHitLink", 1, "extension/MutableSenseWireHitSimTrackerHitLink.h", 40,
                  typeid(::extension::MutableSenseWireHitSimTrackerHitLink), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMutableSenseWireHitSimTrackerHitLink_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MutableSenseWireHitSimTrackerHitLink), alignof(::extension::MutableSenseWireHitSimTrackerHitLink) );
      instance.SetNew(&new_extensioncLcLMutableSenseWireHitSimTrackerHitLink);
      instance.SetNewArray(&newArray_extensioncLcLMutableSenseWireHitSimTrackerHitLink);
      instance.SetDelete(&delete_extensioncLcLMutableSenseWireHitSimTrackerHitLink);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMutableSenseWireHitSimTrackerHitLink);
      instance.SetDestructor(&destruct_extensioncLcLMutableSenseWireHitSimTrackerHitLink);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MutableSenseWireHitSimTrackerHitLink*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MutableSenseWireHitSimTrackerHitLink*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MutableSenseWireHitSimTrackerHitLink*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMutableSenseWireHitSimTrackerHitLink_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MutableSenseWireHitSimTrackerHitLink*>(nullptr))->GetClass();
      extensioncLcLMutableSenseWireHitSimTrackerHitLink_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMutableSenseWireHitSimTrackerHitLink_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLSenseWireHitSimTrackerHitLinkCollection_Dictionary();
   static void extensioncLcLSenseWireHitSimTrackerHitLinkCollection_TClassManip(TClass*);
   static void *new_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(void *p = nullptr);
   static void *newArray_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(Long_t size, void *p);
   static void delete_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(void *p);
   static void deleteArray_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(void *p);
   static void destruct_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::SenseWireHitSimTrackerHitLinkCollection*)
   {
      ::extension::SenseWireHitSimTrackerHitLinkCollection *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::SenseWireHitSimTrackerHitLinkCollection));
      static ::ROOT::TGenericClassInfo 
         instance("extension::SenseWireHitSimTrackerHitLinkCollection", 1, "", 2559,
                  typeid(::extension::SenseWireHitSimTrackerHitLinkCollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLSenseWireHitSimTrackerHitLinkCollection_Dictionary, isa_proxy, 12,
                  sizeof(::extension::SenseWireHitSimTrackerHitLinkCollection), alignof(::extension::SenseWireHitSimTrackerHitLinkCollection) );
      instance.SetNew(&new_extensioncLcLSenseWireHitSimTrackerHitLinkCollection);
      instance.SetNewArray(&newArray_extensioncLcLSenseWireHitSimTrackerHitLinkCollection);
      instance.SetDelete(&delete_extensioncLcLSenseWireHitSimTrackerHitLinkCollection);
      instance.SetDeleteArray(&deleteArray_extensioncLcLSenseWireHitSimTrackerHitLinkCollection);
      instance.SetDestructor(&destruct_extensioncLcLSenseWireHitSimTrackerHitLinkCollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::SenseWireHitSimTrackerHitLinkCollection*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::SenseWireHitSimTrackerHitLinkCollection*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitSimTrackerHitLinkCollection*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLSenseWireHitSimTrackerHitLinkCollection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::SenseWireHitSimTrackerHitLinkCollection*>(nullptr))->GetClass();
      extensioncLcLSenseWireHitSimTrackerHitLinkCollection_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLSenseWireHitSimTrackerHitLinkCollection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLTrackData_Dictionary();
   static void extensioncLcLTrackData_TClassManip(TClass*);
   static void *new_extensioncLcLTrackData(void *p = nullptr);
   static void *newArray_extensioncLcLTrackData(Long_t size, void *p);
   static void delete_extensioncLcLTrackData(void *p);
   static void deleteArray_extensioncLcLTrackData(void *p);
   static void destruct_extensioncLcLTrackData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::TrackData*)
   {
      ::extension::TrackData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::TrackData));
      static ::ROOT::TGenericClassInfo 
         instance("extension::TrackData", 1, "extension/TrackData.h", 19,
                  typeid(::extension::TrackData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLTrackData_Dictionary, isa_proxy, 12,
                  sizeof(::extension::TrackData), alignof(::extension::TrackData) );
      instance.SetNew(&new_extensioncLcLTrackData);
      instance.SetNewArray(&newArray_extensioncLcLTrackData);
      instance.SetDelete(&delete_extensioncLcLTrackData);
      instance.SetDeleteArray(&deleteArray_extensioncLcLTrackData);
      instance.SetDestructor(&destruct_extensioncLcLTrackData);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("extension::TrackData","extension::v1::TrackData"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::TrackData*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::TrackData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::TrackData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLTrackData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::TrackData*>(nullptr))->GetClass();
      extensioncLcLTrackData_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLTrackData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLTrack_Dictionary();
   static void extensioncLcLTrack_TClassManip(TClass*);
   static void *new_extensioncLcLTrack(void *p = nullptr);
   static void *newArray_extensioncLcLTrack(Long_t size, void *p);
   static void delete_extensioncLcLTrack(void *p);
   static void deleteArray_extensioncLcLTrack(void *p);
   static void destruct_extensioncLcLTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::Track*)
   {
      ::extension::Track *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::Track));
      static ::ROOT::TGenericClassInfo 
         instance("extension::Track", 1, "extension/Track.h", 46,
                  typeid(::extension::Track), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLTrack_Dictionary, isa_proxy, 12,
                  sizeof(::extension::Track), alignof(::extension::Track) );
      instance.SetNew(&new_extensioncLcLTrack);
      instance.SetNewArray(&newArray_extensioncLcLTrack);
      instance.SetDelete(&delete_extensioncLcLTrack);
      instance.SetDeleteArray(&deleteArray_extensioncLcLTrack);
      instance.SetDestructor(&destruct_extensioncLcLTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::Track*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::Track*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::Track*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLTrack_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::Track*>(nullptr))->GetClass();
      extensioncLcLTrack_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLTrack_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLMutableTrack_Dictionary();
   static void extensioncLcLMutableTrack_TClassManip(TClass*);
   static void *new_extensioncLcLMutableTrack(void *p = nullptr);
   static void *newArray_extensioncLcLMutableTrack(Long_t size, void *p);
   static void delete_extensioncLcLMutableTrack(void *p);
   static void deleteArray_extensioncLcLMutableTrack(void *p);
   static void destruct_extensioncLcLMutableTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::MutableTrack*)
   {
      ::extension::MutableTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::MutableTrack));
      static ::ROOT::TGenericClassInfo 
         instance("extension::MutableTrack", 1, "extension/MutableTrack.h", 38,
                  typeid(::extension::MutableTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLMutableTrack_Dictionary, isa_proxy, 12,
                  sizeof(::extension::MutableTrack), alignof(::extension::MutableTrack) );
      instance.SetNew(&new_extensioncLcLMutableTrack);
      instance.SetNewArray(&newArray_extensioncLcLMutableTrack);
      instance.SetDelete(&delete_extensioncLcLMutableTrack);
      instance.SetDeleteArray(&deleteArray_extensioncLcLMutableTrack);
      instance.SetDestructor(&destruct_extensioncLcLMutableTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::MutableTrack*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::MutableTrack*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::MutableTrack*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLMutableTrack_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::MutableTrack*>(nullptr))->GetClass();
      extensioncLcLMutableTrack_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLMutableTrack_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *extensioncLcLTrackCollection_Dictionary();
   static void extensioncLcLTrackCollection_TClassManip(TClass*);
   static void *new_extensioncLcLTrackCollection(void *p = nullptr);
   static void *newArray_extensioncLcLTrackCollection(Long_t size, void *p);
   static void delete_extensioncLcLTrackCollection(void *p);
   static void deleteArray_extensioncLcLTrackCollection(void *p);
   static void destruct_extensioncLcLTrackCollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::extension::TrackCollection*)
   {
      ::extension::TrackCollection *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::extension::TrackCollection));
      static ::ROOT::TGenericClassInfo 
         instance("extension::TrackCollection", 1, "", 2954,
                  typeid(::extension::TrackCollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &extensioncLcLTrackCollection_Dictionary, isa_proxy, 12,
                  sizeof(::extension::TrackCollection), alignof(::extension::TrackCollection) );
      instance.SetNew(&new_extensioncLcLTrackCollection);
      instance.SetNewArray(&newArray_extensioncLcLTrackCollection);
      instance.SetDelete(&delete_extensioncLcLTrackCollection);
      instance.SetDeleteArray(&deleteArray_extensioncLcLTrackCollection);
      instance.SetDestructor(&destruct_extensioncLcLTrackCollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::extension::TrackCollection*)
   {
      return GenerateInitInstanceLocal(static_cast<::extension::TrackCollection*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::extension::TrackCollection*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *extensioncLcLTrackCollection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::extension::TrackCollection*>(nullptr))->GetClass();
      extensioncLcLTrackCollection_TClassManip(theClass);
   return theClass;
   }

   static void extensioncLcLTrackCollection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLDriftChamberDigiData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiData : new ::extension::DriftChamberDigiData;
   }
   static void *newArray_extensioncLcLDriftChamberDigiData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiData[nElements] : new ::extension::DriftChamberDigiData[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLDriftChamberDigiData(void *p) {
      delete (static_cast<::extension::DriftChamberDigiData*>(p));
   }
   static void deleteArray_extensioncLcLDriftChamberDigiData(void *p) {
      delete [] (static_cast<::extension::DriftChamberDigiData*>(p));
   }
   static void destruct_extensioncLcLDriftChamberDigiData(void *p) {
      typedef ::extension::DriftChamberDigiData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::DriftChamberDigiData

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLDriftChamberDigi(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigi : new ::extension::DriftChamberDigi;
   }
   static void *newArray_extensioncLcLDriftChamberDigi(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigi[nElements] : new ::extension::DriftChamberDigi[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLDriftChamberDigi(void *p) {
      delete (static_cast<::extension::DriftChamberDigi*>(p));
   }
   static void deleteArray_extensioncLcLDriftChamberDigi(void *p) {
      delete [] (static_cast<::extension::DriftChamberDigi*>(p));
   }
   static void destruct_extensioncLcLDriftChamberDigi(void *p) {
      typedef ::extension::DriftChamberDigi current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::DriftChamberDigi

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMutableDriftChamberDigi(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableDriftChamberDigi : new ::extension::MutableDriftChamberDigi;
   }
   static void *newArray_extensioncLcLMutableDriftChamberDigi(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableDriftChamberDigi[nElements] : new ::extension::MutableDriftChamberDigi[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMutableDriftChamberDigi(void *p) {
      delete (static_cast<::extension::MutableDriftChamberDigi*>(p));
   }
   static void deleteArray_extensioncLcLMutableDriftChamberDigi(void *p) {
      delete [] (static_cast<::extension::MutableDriftChamberDigi*>(p));
   }
   static void destruct_extensioncLcLMutableDriftChamberDigi(void *p) {
      typedef ::extension::MutableDriftChamberDigi current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MutableDriftChamberDigi

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLDriftChamberDigiCollection(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiCollection : new ::extension::DriftChamberDigiCollection;
   }
   static void *newArray_extensioncLcLDriftChamberDigiCollection(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiCollection[nElements] : new ::extension::DriftChamberDigiCollection[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLDriftChamberDigiCollection(void *p) {
      delete (static_cast<::extension::DriftChamberDigiCollection*>(p));
   }
   static void deleteArray_extensioncLcLDriftChamberDigiCollection(void *p) {
      delete [] (static_cast<::extension::DriftChamberDigiCollection*>(p));
   }
   static void destruct_extensioncLcLDriftChamberDigiCollection(void *p) {
      typedef ::extension::DriftChamberDigiCollection current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::DriftChamberDigiCollection

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLDriftChamberDigiLocalData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiLocalData : new ::extension::DriftChamberDigiLocalData;
   }
   static void *newArray_extensioncLcLDriftChamberDigiLocalData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiLocalData[nElements] : new ::extension::DriftChamberDigiLocalData[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLDriftChamberDigiLocalData(void *p) {
      delete (static_cast<::extension::DriftChamberDigiLocalData*>(p));
   }
   static void deleteArray_extensioncLcLDriftChamberDigiLocalData(void *p) {
      delete [] (static_cast<::extension::DriftChamberDigiLocalData*>(p));
   }
   static void destruct_extensioncLcLDriftChamberDigiLocalData(void *p) {
      typedef ::extension::DriftChamberDigiLocalData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::DriftChamberDigiLocalData

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLDriftChamberDigiLocal(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiLocal : new ::extension::DriftChamberDigiLocal;
   }
   static void *newArray_extensioncLcLDriftChamberDigiLocal(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiLocal[nElements] : new ::extension::DriftChamberDigiLocal[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLDriftChamberDigiLocal(void *p) {
      delete (static_cast<::extension::DriftChamberDigiLocal*>(p));
   }
   static void deleteArray_extensioncLcLDriftChamberDigiLocal(void *p) {
      delete [] (static_cast<::extension::DriftChamberDigiLocal*>(p));
   }
   static void destruct_extensioncLcLDriftChamberDigiLocal(void *p) {
      typedef ::extension::DriftChamberDigiLocal current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::DriftChamberDigiLocal

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMutableDriftChamberDigiLocal(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableDriftChamberDigiLocal : new ::extension::MutableDriftChamberDigiLocal;
   }
   static void *newArray_extensioncLcLMutableDriftChamberDigiLocal(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableDriftChamberDigiLocal[nElements] : new ::extension::MutableDriftChamberDigiLocal[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMutableDriftChamberDigiLocal(void *p) {
      delete (static_cast<::extension::MutableDriftChamberDigiLocal*>(p));
   }
   static void deleteArray_extensioncLcLMutableDriftChamberDigiLocal(void *p) {
      delete [] (static_cast<::extension::MutableDriftChamberDigiLocal*>(p));
   }
   static void destruct_extensioncLcLMutableDriftChamberDigiLocal(void *p) {
      typedef ::extension::MutableDriftChamberDigiLocal current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MutableDriftChamberDigiLocal

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLDriftChamberDigiLocalCollection(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiLocalCollection : new ::extension::DriftChamberDigiLocalCollection;
   }
   static void *newArray_extensioncLcLDriftChamberDigiLocalCollection(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiLocalCollection[nElements] : new ::extension::DriftChamberDigiLocalCollection[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLDriftChamberDigiLocalCollection(void *p) {
      delete (static_cast<::extension::DriftChamberDigiLocalCollection*>(p));
   }
   static void deleteArray_extensioncLcLDriftChamberDigiLocalCollection(void *p) {
      delete [] (static_cast<::extension::DriftChamberDigiLocalCollection*>(p));
   }
   static void destruct_extensioncLcLDriftChamberDigiLocalCollection(void *p) {
      typedef ::extension::DriftChamberDigiLocalCollection current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::DriftChamberDigiLocalCollection

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLDriftChamberDigiV2Data(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiV2Data : new ::extension::DriftChamberDigiV2Data;
   }
   static void *newArray_extensioncLcLDriftChamberDigiV2Data(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiV2Data[nElements] : new ::extension::DriftChamberDigiV2Data[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLDriftChamberDigiV2Data(void *p) {
      delete (static_cast<::extension::DriftChamberDigiV2Data*>(p));
   }
   static void deleteArray_extensioncLcLDriftChamberDigiV2Data(void *p) {
      delete [] (static_cast<::extension::DriftChamberDigiV2Data*>(p));
   }
   static void destruct_extensioncLcLDriftChamberDigiV2Data(void *p) {
      typedef ::extension::DriftChamberDigiV2Data current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::DriftChamberDigiV2Data

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLDriftChamberDigiV2(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiV2 : new ::extension::DriftChamberDigiV2;
   }
   static void *newArray_extensioncLcLDriftChamberDigiV2(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiV2[nElements] : new ::extension::DriftChamberDigiV2[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLDriftChamberDigiV2(void *p) {
      delete (static_cast<::extension::DriftChamberDigiV2*>(p));
   }
   static void deleteArray_extensioncLcLDriftChamberDigiV2(void *p) {
      delete [] (static_cast<::extension::DriftChamberDigiV2*>(p));
   }
   static void destruct_extensioncLcLDriftChamberDigiV2(void *p) {
      typedef ::extension::DriftChamberDigiV2 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::DriftChamberDigiV2

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMutableDriftChamberDigiV2(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableDriftChamberDigiV2 : new ::extension::MutableDriftChamberDigiV2;
   }
   static void *newArray_extensioncLcLMutableDriftChamberDigiV2(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableDriftChamberDigiV2[nElements] : new ::extension::MutableDriftChamberDigiV2[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMutableDriftChamberDigiV2(void *p) {
      delete (static_cast<::extension::MutableDriftChamberDigiV2*>(p));
   }
   static void deleteArray_extensioncLcLMutableDriftChamberDigiV2(void *p) {
      delete [] (static_cast<::extension::MutableDriftChamberDigiV2*>(p));
   }
   static void destruct_extensioncLcLMutableDriftChamberDigiV2(void *p) {
      typedef ::extension::MutableDriftChamberDigiV2 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MutableDriftChamberDigiV2

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLDriftChamberDigiV2Collection(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiV2Collection : new ::extension::DriftChamberDigiV2Collection;
   }
   static void *newArray_extensioncLcLDriftChamberDigiV2Collection(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::DriftChamberDigiV2Collection[nElements] : new ::extension::DriftChamberDigiV2Collection[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLDriftChamberDigiV2Collection(void *p) {
      delete (static_cast<::extension::DriftChamberDigiV2Collection*>(p));
   }
   static void deleteArray_extensioncLcLDriftChamberDigiV2Collection(void *p) {
      delete [] (static_cast<::extension::DriftChamberDigiV2Collection*>(p));
   }
   static void destruct_extensioncLcLDriftChamberDigiV2Collection(void *p) {
      typedef ::extension::DriftChamberDigiV2Collection current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::DriftChamberDigiV2Collection

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMCRecoDriftChamberDigiAssociationData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiAssociationData : new ::extension::MCRecoDriftChamberDigiAssociationData;
   }
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiAssociationData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiAssociationData[nElements] : new ::extension::MCRecoDriftChamberDigiAssociationData[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMCRecoDriftChamberDigiAssociationData(void *p) {
      delete (static_cast<::extension::MCRecoDriftChamberDigiAssociationData*>(p));
   }
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiAssociationData(void *p) {
      delete [] (static_cast<::extension::MCRecoDriftChamberDigiAssociationData*>(p));
   }
   static void destruct_extensioncLcLMCRecoDriftChamberDigiAssociationData(void *p) {
      typedef ::extension::MCRecoDriftChamberDigiAssociationData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MCRecoDriftChamberDigiAssociationData

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMCRecoDriftChamberDigiAssociation(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiAssociation : new ::extension::MCRecoDriftChamberDigiAssociation;
   }
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiAssociation(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiAssociation[nElements] : new ::extension::MCRecoDriftChamberDigiAssociation[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMCRecoDriftChamberDigiAssociation(void *p) {
      delete (static_cast<::extension::MCRecoDriftChamberDigiAssociation*>(p));
   }
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiAssociation(void *p) {
      delete [] (static_cast<::extension::MCRecoDriftChamberDigiAssociation*>(p));
   }
   static void destruct_extensioncLcLMCRecoDriftChamberDigiAssociation(void *p) {
      typedef ::extension::MCRecoDriftChamberDigiAssociation current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MCRecoDriftChamberDigiAssociation

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableMCRecoDriftChamberDigiAssociation : new ::extension::MutableMCRecoDriftChamberDigiAssociation;
   }
   static void *newArray_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableMCRecoDriftChamberDigiAssociation[nElements] : new ::extension::MutableMCRecoDriftChamberDigiAssociation[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(void *p) {
      delete (static_cast<::extension::MutableMCRecoDriftChamberDigiAssociation*>(p));
   }
   static void deleteArray_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(void *p) {
      delete [] (static_cast<::extension::MutableMCRecoDriftChamberDigiAssociation*>(p));
   }
   static void destruct_extensioncLcLMutableMCRecoDriftChamberDigiAssociation(void *p) {
      typedef ::extension::MutableMCRecoDriftChamberDigiAssociation current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MutableMCRecoDriftChamberDigiAssociation

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiAssociationCollection : new ::extension::MCRecoDriftChamberDigiAssociationCollection;
   }
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiAssociationCollection[nElements] : new ::extension::MCRecoDriftChamberDigiAssociationCollection[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(void *p) {
      delete (static_cast<::extension::MCRecoDriftChamberDigiAssociationCollection*>(p));
   }
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(void *p) {
      delete [] (static_cast<::extension::MCRecoDriftChamberDigiAssociationCollection*>(p));
   }
   static void destruct_extensioncLcLMCRecoDriftChamberDigiAssociationCollection(void *p) {
      typedef ::extension::MCRecoDriftChamberDigiAssociationCollection current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MCRecoDriftChamberDigiAssociationCollection

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiV2AssociationData : new ::extension::MCRecoDriftChamberDigiV2AssociationData;
   }
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiV2AssociationData[nElements] : new ::extension::MCRecoDriftChamberDigiV2AssociationData[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(void *p) {
      delete (static_cast<::extension::MCRecoDriftChamberDigiV2AssociationData*>(p));
   }
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(void *p) {
      delete [] (static_cast<::extension::MCRecoDriftChamberDigiV2AssociationData*>(p));
   }
   static void destruct_extensioncLcLMCRecoDriftChamberDigiV2AssociationData(void *p) {
      typedef ::extension::MCRecoDriftChamberDigiV2AssociationData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MCRecoDriftChamberDigiV2AssociationData

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMCRecoDriftChamberDigiV2Association(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiV2Association : new ::extension::MCRecoDriftChamberDigiV2Association;
   }
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiV2Association(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiV2Association[nElements] : new ::extension::MCRecoDriftChamberDigiV2Association[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMCRecoDriftChamberDigiV2Association(void *p) {
      delete (static_cast<::extension::MCRecoDriftChamberDigiV2Association*>(p));
   }
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiV2Association(void *p) {
      delete [] (static_cast<::extension::MCRecoDriftChamberDigiV2Association*>(p));
   }
   static void destruct_extensioncLcLMCRecoDriftChamberDigiV2Association(void *p) {
      typedef ::extension::MCRecoDriftChamberDigiV2Association current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MCRecoDriftChamberDigiV2Association

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableMCRecoDriftChamberDigiV2Association : new ::extension::MutableMCRecoDriftChamberDigiV2Association;
   }
   static void *newArray_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableMCRecoDriftChamberDigiV2Association[nElements] : new ::extension::MutableMCRecoDriftChamberDigiV2Association[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(void *p) {
      delete (static_cast<::extension::MutableMCRecoDriftChamberDigiV2Association*>(p));
   }
   static void deleteArray_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(void *p) {
      delete [] (static_cast<::extension::MutableMCRecoDriftChamberDigiV2Association*>(p));
   }
   static void destruct_extensioncLcLMutableMCRecoDriftChamberDigiV2Association(void *p) {
      typedef ::extension::MutableMCRecoDriftChamberDigiV2Association current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MutableMCRecoDriftChamberDigiV2Association

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiV2AssociationCollection : new ::extension::MCRecoDriftChamberDigiV2AssociationCollection;
   }
   static void *newArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MCRecoDriftChamberDigiV2AssociationCollection[nElements] : new ::extension::MCRecoDriftChamberDigiV2AssociationCollection[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(void *p) {
      delete (static_cast<::extension::MCRecoDriftChamberDigiV2AssociationCollection*>(p));
   }
   static void deleteArray_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(void *p) {
      delete [] (static_cast<::extension::MCRecoDriftChamberDigiV2AssociationCollection*>(p));
   }
   static void destruct_extensioncLcLMCRecoDriftChamberDigiV2AssociationCollection(void *p) {
      typedef ::extension::MCRecoDriftChamberDigiV2AssociationCollection current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MCRecoDriftChamberDigiV2AssociationCollection

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLSenseWireHitData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitData : new ::extension::SenseWireHitData;
   }
   static void *newArray_extensioncLcLSenseWireHitData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitData[nElements] : new ::extension::SenseWireHitData[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLSenseWireHitData(void *p) {
      delete (static_cast<::extension::SenseWireHitData*>(p));
   }
   static void deleteArray_extensioncLcLSenseWireHitData(void *p) {
      delete [] (static_cast<::extension::SenseWireHitData*>(p));
   }
   static void destruct_extensioncLcLSenseWireHitData(void *p) {
      typedef ::extension::SenseWireHitData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::SenseWireHitData

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLSenseWireHit(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHit : new ::extension::SenseWireHit;
   }
   static void *newArray_extensioncLcLSenseWireHit(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHit[nElements] : new ::extension::SenseWireHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLSenseWireHit(void *p) {
      delete (static_cast<::extension::SenseWireHit*>(p));
   }
   static void deleteArray_extensioncLcLSenseWireHit(void *p) {
      delete [] (static_cast<::extension::SenseWireHit*>(p));
   }
   static void destruct_extensioncLcLSenseWireHit(void *p) {
      typedef ::extension::SenseWireHit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::SenseWireHit

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMutableSenseWireHit(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableSenseWireHit : new ::extension::MutableSenseWireHit;
   }
   static void *newArray_extensioncLcLMutableSenseWireHit(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableSenseWireHit[nElements] : new ::extension::MutableSenseWireHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMutableSenseWireHit(void *p) {
      delete (static_cast<::extension::MutableSenseWireHit*>(p));
   }
   static void deleteArray_extensioncLcLMutableSenseWireHit(void *p) {
      delete [] (static_cast<::extension::MutableSenseWireHit*>(p));
   }
   static void destruct_extensioncLcLMutableSenseWireHit(void *p) {
      typedef ::extension::MutableSenseWireHit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MutableSenseWireHit

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLSenseWireHitCollection(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitCollection : new ::extension::SenseWireHitCollection;
   }
   static void *newArray_extensioncLcLSenseWireHitCollection(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitCollection[nElements] : new ::extension::SenseWireHitCollection[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLSenseWireHitCollection(void *p) {
      delete (static_cast<::extension::SenseWireHitCollection*>(p));
   }
   static void deleteArray_extensioncLcLSenseWireHitCollection(void *p) {
      delete [] (static_cast<::extension::SenseWireHitCollection*>(p));
   }
   static void destruct_extensioncLcLSenseWireHitCollection(void *p) {
      typedef ::extension::SenseWireHitCollection current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::SenseWireHitCollection

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLSenseWireHitSimTrackerHitLinkData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitSimTrackerHitLinkData : new ::extension::SenseWireHitSimTrackerHitLinkData;
   }
   static void *newArray_extensioncLcLSenseWireHitSimTrackerHitLinkData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitSimTrackerHitLinkData[nElements] : new ::extension::SenseWireHitSimTrackerHitLinkData[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLSenseWireHitSimTrackerHitLinkData(void *p) {
      delete (static_cast<::extension::SenseWireHitSimTrackerHitLinkData*>(p));
   }
   static void deleteArray_extensioncLcLSenseWireHitSimTrackerHitLinkData(void *p) {
      delete [] (static_cast<::extension::SenseWireHitSimTrackerHitLinkData*>(p));
   }
   static void destruct_extensioncLcLSenseWireHitSimTrackerHitLinkData(void *p) {
      typedef ::extension::SenseWireHitSimTrackerHitLinkData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::SenseWireHitSimTrackerHitLinkData

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLSenseWireHitSimTrackerHitLink(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitSimTrackerHitLink : new ::extension::SenseWireHitSimTrackerHitLink;
   }
   static void *newArray_extensioncLcLSenseWireHitSimTrackerHitLink(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitSimTrackerHitLink[nElements] : new ::extension::SenseWireHitSimTrackerHitLink[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLSenseWireHitSimTrackerHitLink(void *p) {
      delete (static_cast<::extension::SenseWireHitSimTrackerHitLink*>(p));
   }
   static void deleteArray_extensioncLcLSenseWireHitSimTrackerHitLink(void *p) {
      delete [] (static_cast<::extension::SenseWireHitSimTrackerHitLink*>(p));
   }
   static void destruct_extensioncLcLSenseWireHitSimTrackerHitLink(void *p) {
      typedef ::extension::SenseWireHitSimTrackerHitLink current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::SenseWireHitSimTrackerHitLink

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMutableSenseWireHitSimTrackerHitLink(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableSenseWireHitSimTrackerHitLink : new ::extension::MutableSenseWireHitSimTrackerHitLink;
   }
   static void *newArray_extensioncLcLMutableSenseWireHitSimTrackerHitLink(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableSenseWireHitSimTrackerHitLink[nElements] : new ::extension::MutableSenseWireHitSimTrackerHitLink[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMutableSenseWireHitSimTrackerHitLink(void *p) {
      delete (static_cast<::extension::MutableSenseWireHitSimTrackerHitLink*>(p));
   }
   static void deleteArray_extensioncLcLMutableSenseWireHitSimTrackerHitLink(void *p) {
      delete [] (static_cast<::extension::MutableSenseWireHitSimTrackerHitLink*>(p));
   }
   static void destruct_extensioncLcLMutableSenseWireHitSimTrackerHitLink(void *p) {
      typedef ::extension::MutableSenseWireHitSimTrackerHitLink current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MutableSenseWireHitSimTrackerHitLink

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitSimTrackerHitLinkCollection : new ::extension::SenseWireHitSimTrackerHitLinkCollection;
   }
   static void *newArray_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::SenseWireHitSimTrackerHitLinkCollection[nElements] : new ::extension::SenseWireHitSimTrackerHitLinkCollection[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(void *p) {
      delete (static_cast<::extension::SenseWireHitSimTrackerHitLinkCollection*>(p));
   }
   static void deleteArray_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(void *p) {
      delete [] (static_cast<::extension::SenseWireHitSimTrackerHitLinkCollection*>(p));
   }
   static void destruct_extensioncLcLSenseWireHitSimTrackerHitLinkCollection(void *p) {
      typedef ::extension::SenseWireHitSimTrackerHitLinkCollection current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::SenseWireHitSimTrackerHitLinkCollection

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLTrackData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::TrackData : new ::extension::TrackData;
   }
   static void *newArray_extensioncLcLTrackData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::TrackData[nElements] : new ::extension::TrackData[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLTrackData(void *p) {
      delete (static_cast<::extension::TrackData*>(p));
   }
   static void deleteArray_extensioncLcLTrackData(void *p) {
      delete [] (static_cast<::extension::TrackData*>(p));
   }
   static void destruct_extensioncLcLTrackData(void *p) {
      typedef ::extension::TrackData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::TrackData

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLTrack(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::Track : new ::extension::Track;
   }
   static void *newArray_extensioncLcLTrack(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::Track[nElements] : new ::extension::Track[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLTrack(void *p) {
      delete (static_cast<::extension::Track*>(p));
   }
   static void deleteArray_extensioncLcLTrack(void *p) {
      delete [] (static_cast<::extension::Track*>(p));
   }
   static void destruct_extensioncLcLTrack(void *p) {
      typedef ::extension::Track current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::Track

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLMutableTrack(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableTrack : new ::extension::MutableTrack;
   }
   static void *newArray_extensioncLcLMutableTrack(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::MutableTrack[nElements] : new ::extension::MutableTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLMutableTrack(void *p) {
      delete (static_cast<::extension::MutableTrack*>(p));
   }
   static void deleteArray_extensioncLcLMutableTrack(void *p) {
      delete [] (static_cast<::extension::MutableTrack*>(p));
   }
   static void destruct_extensioncLcLMutableTrack(void *p) {
      typedef ::extension::MutableTrack current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::MutableTrack

namespace ROOT {
   // Wrappers around operator new
   static void *new_extensioncLcLTrackCollection(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::TrackCollection : new ::extension::TrackCollection;
   }
   static void *newArray_extensioncLcLTrackCollection(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::extension::TrackCollection[nElements] : new ::extension::TrackCollection[nElements];
   }
   // Wrapper around operator delete
   static void delete_extensioncLcLTrackCollection(void *p) {
      delete (static_cast<::extension::TrackCollection*>(p));
   }
   static void deleteArray_extensioncLcLTrackCollection(void *p) {
      delete [] (static_cast<::extension::TrackCollection*>(p));
   }
   static void destruct_extensioncLcLTrackCollection(void *p) {
      typedef ::extension::TrackCollection current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::extension::TrackCollection

namespace ROOT {
   static TClass *vectorlEextensioncLcLTrackDatagR_Dictionary();
   static void vectorlEextensioncLcLTrackDatagR_TClassManip(TClass*);
   static void *new_vectorlEextensioncLcLTrackDatagR(void *p = nullptr);
   static void *newArray_vectorlEextensioncLcLTrackDatagR(Long_t size, void *p);
   static void delete_vectorlEextensioncLcLTrackDatagR(void *p);
   static void deleteArray_vectorlEextensioncLcLTrackDatagR(void *p);
   static void destruct_vectorlEextensioncLcLTrackDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<extension::TrackData>*)
   {
      vector<extension::TrackData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<extension::TrackData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<extension::TrackData>", -2, "vector", 460,
                  typeid(vector<extension::TrackData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEextensioncLcLTrackDatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<extension::TrackData>), alignof(vector<extension::TrackData>) );
      instance.SetNew(&new_vectorlEextensioncLcLTrackDatagR);
      instance.SetNewArray(&newArray_vectorlEextensioncLcLTrackDatagR);
      instance.SetDelete(&delete_vectorlEextensioncLcLTrackDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEextensioncLcLTrackDatagR);
      instance.SetDestructor(&destruct_vectorlEextensioncLcLTrackDatagR);
      static_assert(alignof(vector<extension::TrackData>::value_type) <= 4096,
          "Class with alignment strictly greater than 4096 are currently not supported in CollectionProxy. "
          "Please report this case to the developers");
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<extension::TrackData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<extension::TrackData>","std::vector<extension::TrackData, std::allocator<extension::TrackData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<extension::TrackData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEextensioncLcLTrackDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<extension::TrackData>*>(nullptr))->GetClass();
      vectorlEextensioncLcLTrackDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEextensioncLcLTrackDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEextensioncLcLTrackDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::TrackData> : new vector<extension::TrackData>;
   }
   static void *newArray_vectorlEextensioncLcLTrackDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::TrackData>[nElements] : new vector<extension::TrackData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEextensioncLcLTrackDatagR(void *p) {
      delete (static_cast<vector<extension::TrackData>*>(p));
   }
   static void deleteArray_vectorlEextensioncLcLTrackDatagR(void *p) {
      delete [] (static_cast<vector<extension::TrackData>*>(p));
   }
   static void destruct_vectorlEextensioncLcLTrackDatagR(void *p) {
      typedef vector<extension::TrackData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<extension::TrackData>

namespace ROOT {
   static TClass *vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR_Dictionary();
   static void vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR_TClassManip(TClass*);
   static void *new_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(void *p = nullptr);
   static void *newArray_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(Long_t size, void *p);
   static void delete_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(void *p);
   static void deleteArray_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(void *p);
   static void destruct_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<extension::SenseWireHitSimTrackerHitLinkData>*)
   {
      vector<extension::SenseWireHitSimTrackerHitLinkData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<extension::SenseWireHitSimTrackerHitLinkData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<extension::SenseWireHitSimTrackerHitLinkData>", -2, "vector", 460,
                  typeid(vector<extension::SenseWireHitSimTrackerHitLinkData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<extension::SenseWireHitSimTrackerHitLinkData>), alignof(vector<extension::SenseWireHitSimTrackerHitLinkData>) );
      instance.SetNew(&new_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR);
      instance.SetNewArray(&newArray_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR);
      instance.SetDelete(&delete_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR);
      instance.SetDestructor(&destruct_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR);
      static_assert(alignof(vector<extension::SenseWireHitSimTrackerHitLinkData>::value_type) <= 4096,
          "Class with alignment strictly greater than 4096 are currently not supported in CollectionProxy. "
          "Please report this case to the developers");
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<extension::SenseWireHitSimTrackerHitLinkData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<extension::SenseWireHitSimTrackerHitLinkData>","std::vector<extension::SenseWireHitSimTrackerHitLinkData, std::allocator<extension::SenseWireHitSimTrackerHitLinkData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<extension::SenseWireHitSimTrackerHitLinkData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<extension::SenseWireHitSimTrackerHitLinkData>*>(nullptr))->GetClass();
      vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::SenseWireHitSimTrackerHitLinkData> : new vector<extension::SenseWireHitSimTrackerHitLinkData>;
   }
   static void *newArray_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::SenseWireHitSimTrackerHitLinkData>[nElements] : new vector<extension::SenseWireHitSimTrackerHitLinkData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(void *p) {
      delete (static_cast<vector<extension::SenseWireHitSimTrackerHitLinkData>*>(p));
   }
   static void deleteArray_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(void *p) {
      delete [] (static_cast<vector<extension::SenseWireHitSimTrackerHitLinkData>*>(p));
   }
   static void destruct_vectorlEextensioncLcLSenseWireHitSimTrackerHitLinkDatagR(void *p) {
      typedef vector<extension::SenseWireHitSimTrackerHitLinkData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<extension::SenseWireHitSimTrackerHitLinkData>

namespace ROOT {
   static TClass *vectorlEextensioncLcLSenseWireHitDatagR_Dictionary();
   static void vectorlEextensioncLcLSenseWireHitDatagR_TClassManip(TClass*);
   static void *new_vectorlEextensioncLcLSenseWireHitDatagR(void *p = nullptr);
   static void *newArray_vectorlEextensioncLcLSenseWireHitDatagR(Long_t size, void *p);
   static void delete_vectorlEextensioncLcLSenseWireHitDatagR(void *p);
   static void deleteArray_vectorlEextensioncLcLSenseWireHitDatagR(void *p);
   static void destruct_vectorlEextensioncLcLSenseWireHitDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<extension::SenseWireHitData>*)
   {
      vector<extension::SenseWireHitData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<extension::SenseWireHitData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<extension::SenseWireHitData>", -2, "vector", 460,
                  typeid(vector<extension::SenseWireHitData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEextensioncLcLSenseWireHitDatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<extension::SenseWireHitData>), alignof(vector<extension::SenseWireHitData>) );
      instance.SetNew(&new_vectorlEextensioncLcLSenseWireHitDatagR);
      instance.SetNewArray(&newArray_vectorlEextensioncLcLSenseWireHitDatagR);
      instance.SetDelete(&delete_vectorlEextensioncLcLSenseWireHitDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEextensioncLcLSenseWireHitDatagR);
      instance.SetDestructor(&destruct_vectorlEextensioncLcLSenseWireHitDatagR);
      static_assert(alignof(vector<extension::SenseWireHitData>::value_type) <= 4096,
          "Class with alignment strictly greater than 4096 are currently not supported in CollectionProxy. "
          "Please report this case to the developers");
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<extension::SenseWireHitData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<extension::SenseWireHitData>","std::vector<extension::SenseWireHitData, std::allocator<extension::SenseWireHitData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<extension::SenseWireHitData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEextensioncLcLSenseWireHitDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<extension::SenseWireHitData>*>(nullptr))->GetClass();
      vectorlEextensioncLcLSenseWireHitDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEextensioncLcLSenseWireHitDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEextensioncLcLSenseWireHitDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::SenseWireHitData> : new vector<extension::SenseWireHitData>;
   }
   static void *newArray_vectorlEextensioncLcLSenseWireHitDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::SenseWireHitData>[nElements] : new vector<extension::SenseWireHitData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEextensioncLcLSenseWireHitDatagR(void *p) {
      delete (static_cast<vector<extension::SenseWireHitData>*>(p));
   }
   static void deleteArray_vectorlEextensioncLcLSenseWireHitDatagR(void *p) {
      delete [] (static_cast<vector<extension::SenseWireHitData>*>(p));
   }
   static void destruct_vectorlEextensioncLcLSenseWireHitDatagR(void *p) {
      typedef vector<extension::SenseWireHitData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<extension::SenseWireHitData>

namespace ROOT {
   static TClass *vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR_Dictionary();
   static void vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR_TClassManip(TClass*);
   static void *new_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(void *p = nullptr);
   static void *newArray_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(Long_t size, void *p);
   static void delete_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(void *p);
   static void deleteArray_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(void *p);
   static void destruct_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<extension::MCRecoDriftChamberDigiV2AssociationData>*)
   {
      vector<extension::MCRecoDriftChamberDigiV2AssociationData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<extension::MCRecoDriftChamberDigiV2AssociationData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<extension::MCRecoDriftChamberDigiV2AssociationData>", -2, "vector", 460,
                  typeid(vector<extension::MCRecoDriftChamberDigiV2AssociationData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<extension::MCRecoDriftChamberDigiV2AssociationData>), alignof(vector<extension::MCRecoDriftChamberDigiV2AssociationData>) );
      instance.SetNew(&new_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR);
      instance.SetNewArray(&newArray_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR);
      instance.SetDelete(&delete_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR);
      instance.SetDestructor(&destruct_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR);
      static_assert(alignof(vector<extension::MCRecoDriftChamberDigiV2AssociationData>::value_type) <= 4096,
          "Class with alignment strictly greater than 4096 are currently not supported in CollectionProxy. "
          "Please report this case to the developers");
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<extension::MCRecoDriftChamberDigiV2AssociationData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<extension::MCRecoDriftChamberDigiV2AssociationData>","std::vector<extension::MCRecoDriftChamberDigiV2AssociationData, std::allocator<extension::MCRecoDriftChamberDigiV2AssociationData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<extension::MCRecoDriftChamberDigiV2AssociationData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<extension::MCRecoDriftChamberDigiV2AssociationData>*>(nullptr))->GetClass();
      vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::MCRecoDriftChamberDigiV2AssociationData> : new vector<extension::MCRecoDriftChamberDigiV2AssociationData>;
   }
   static void *newArray_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::MCRecoDriftChamberDigiV2AssociationData>[nElements] : new vector<extension::MCRecoDriftChamberDigiV2AssociationData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(void *p) {
      delete (static_cast<vector<extension::MCRecoDriftChamberDigiV2AssociationData>*>(p));
   }
   static void deleteArray_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(void *p) {
      delete [] (static_cast<vector<extension::MCRecoDriftChamberDigiV2AssociationData>*>(p));
   }
   static void destruct_vectorlEextensioncLcLMCRecoDriftChamberDigiV2AssociationDatagR(void *p) {
      typedef vector<extension::MCRecoDriftChamberDigiV2AssociationData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<extension::MCRecoDriftChamberDigiV2AssociationData>

namespace ROOT {
   static TClass *vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR_Dictionary();
   static void vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR_TClassManip(TClass*);
   static void *new_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(void *p = nullptr);
   static void *newArray_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(Long_t size, void *p);
   static void delete_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(void *p);
   static void deleteArray_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(void *p);
   static void destruct_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<extension::MCRecoDriftChamberDigiAssociationData>*)
   {
      vector<extension::MCRecoDriftChamberDigiAssociationData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<extension::MCRecoDriftChamberDigiAssociationData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<extension::MCRecoDriftChamberDigiAssociationData>", -2, "vector", 460,
                  typeid(vector<extension::MCRecoDriftChamberDigiAssociationData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<extension::MCRecoDriftChamberDigiAssociationData>), alignof(vector<extension::MCRecoDriftChamberDigiAssociationData>) );
      instance.SetNew(&new_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR);
      instance.SetNewArray(&newArray_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR);
      instance.SetDelete(&delete_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR);
      instance.SetDestructor(&destruct_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR);
      static_assert(alignof(vector<extension::MCRecoDriftChamberDigiAssociationData>::value_type) <= 4096,
          "Class with alignment strictly greater than 4096 are currently not supported in CollectionProxy. "
          "Please report this case to the developers");
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<extension::MCRecoDriftChamberDigiAssociationData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<extension::MCRecoDriftChamberDigiAssociationData>","std::vector<extension::MCRecoDriftChamberDigiAssociationData, std::allocator<extension::MCRecoDriftChamberDigiAssociationData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<extension::MCRecoDriftChamberDigiAssociationData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<extension::MCRecoDriftChamberDigiAssociationData>*>(nullptr))->GetClass();
      vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::MCRecoDriftChamberDigiAssociationData> : new vector<extension::MCRecoDriftChamberDigiAssociationData>;
   }
   static void *newArray_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::MCRecoDriftChamberDigiAssociationData>[nElements] : new vector<extension::MCRecoDriftChamberDigiAssociationData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(void *p) {
      delete (static_cast<vector<extension::MCRecoDriftChamberDigiAssociationData>*>(p));
   }
   static void deleteArray_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(void *p) {
      delete [] (static_cast<vector<extension::MCRecoDriftChamberDigiAssociationData>*>(p));
   }
   static void destruct_vectorlEextensioncLcLMCRecoDriftChamberDigiAssociationDatagR(void *p) {
      typedef vector<extension::MCRecoDriftChamberDigiAssociationData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<extension::MCRecoDriftChamberDigiAssociationData>

namespace ROOT {
   static TClass *vectorlEextensioncLcLDriftChamberDigiV2DatagR_Dictionary();
   static void vectorlEextensioncLcLDriftChamberDigiV2DatagR_TClassManip(TClass*);
   static void *new_vectorlEextensioncLcLDriftChamberDigiV2DatagR(void *p = nullptr);
   static void *newArray_vectorlEextensioncLcLDriftChamberDigiV2DatagR(Long_t size, void *p);
   static void delete_vectorlEextensioncLcLDriftChamberDigiV2DatagR(void *p);
   static void deleteArray_vectorlEextensioncLcLDriftChamberDigiV2DatagR(void *p);
   static void destruct_vectorlEextensioncLcLDriftChamberDigiV2DatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<extension::DriftChamberDigiV2Data>*)
   {
      vector<extension::DriftChamberDigiV2Data> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<extension::DriftChamberDigiV2Data>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<extension::DriftChamberDigiV2Data>", -2, "vector", 460,
                  typeid(vector<extension::DriftChamberDigiV2Data>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEextensioncLcLDriftChamberDigiV2DatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<extension::DriftChamberDigiV2Data>), alignof(vector<extension::DriftChamberDigiV2Data>) );
      instance.SetNew(&new_vectorlEextensioncLcLDriftChamberDigiV2DatagR);
      instance.SetNewArray(&newArray_vectorlEextensioncLcLDriftChamberDigiV2DatagR);
      instance.SetDelete(&delete_vectorlEextensioncLcLDriftChamberDigiV2DatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEextensioncLcLDriftChamberDigiV2DatagR);
      instance.SetDestructor(&destruct_vectorlEextensioncLcLDriftChamberDigiV2DatagR);
      static_assert(alignof(vector<extension::DriftChamberDigiV2Data>::value_type) <= 4096,
          "Class with alignment strictly greater than 4096 are currently not supported in CollectionProxy. "
          "Please report this case to the developers");
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<extension::DriftChamberDigiV2Data> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<extension::DriftChamberDigiV2Data>","std::vector<extension::DriftChamberDigiV2Data, std::allocator<extension::DriftChamberDigiV2Data> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<extension::DriftChamberDigiV2Data>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEextensioncLcLDriftChamberDigiV2DatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<extension::DriftChamberDigiV2Data>*>(nullptr))->GetClass();
      vectorlEextensioncLcLDriftChamberDigiV2DatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEextensioncLcLDriftChamberDigiV2DatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEextensioncLcLDriftChamberDigiV2DatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::DriftChamberDigiV2Data> : new vector<extension::DriftChamberDigiV2Data>;
   }
   static void *newArray_vectorlEextensioncLcLDriftChamberDigiV2DatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::DriftChamberDigiV2Data>[nElements] : new vector<extension::DriftChamberDigiV2Data>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEextensioncLcLDriftChamberDigiV2DatagR(void *p) {
      delete (static_cast<vector<extension::DriftChamberDigiV2Data>*>(p));
   }
   static void deleteArray_vectorlEextensioncLcLDriftChamberDigiV2DatagR(void *p) {
      delete [] (static_cast<vector<extension::DriftChamberDigiV2Data>*>(p));
   }
   static void destruct_vectorlEextensioncLcLDriftChamberDigiV2DatagR(void *p) {
      typedef vector<extension::DriftChamberDigiV2Data> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<extension::DriftChamberDigiV2Data>

namespace ROOT {
   static TClass *vectorlEextensioncLcLDriftChamberDigiLocalDatagR_Dictionary();
   static void vectorlEextensioncLcLDriftChamberDigiLocalDatagR_TClassManip(TClass*);
   static void *new_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(void *p = nullptr);
   static void *newArray_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(Long_t size, void *p);
   static void delete_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(void *p);
   static void deleteArray_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(void *p);
   static void destruct_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<extension::DriftChamberDigiLocalData>*)
   {
      vector<extension::DriftChamberDigiLocalData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<extension::DriftChamberDigiLocalData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<extension::DriftChamberDigiLocalData>", -2, "vector", 460,
                  typeid(vector<extension::DriftChamberDigiLocalData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEextensioncLcLDriftChamberDigiLocalDatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<extension::DriftChamberDigiLocalData>), alignof(vector<extension::DriftChamberDigiLocalData>) );
      instance.SetNew(&new_vectorlEextensioncLcLDriftChamberDigiLocalDatagR);
      instance.SetNewArray(&newArray_vectorlEextensioncLcLDriftChamberDigiLocalDatagR);
      instance.SetDelete(&delete_vectorlEextensioncLcLDriftChamberDigiLocalDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEextensioncLcLDriftChamberDigiLocalDatagR);
      instance.SetDestructor(&destruct_vectorlEextensioncLcLDriftChamberDigiLocalDatagR);
      static_assert(alignof(vector<extension::DriftChamberDigiLocalData>::value_type) <= 4096,
          "Class with alignment strictly greater than 4096 are currently not supported in CollectionProxy. "
          "Please report this case to the developers");
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<extension::DriftChamberDigiLocalData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<extension::DriftChamberDigiLocalData>","std::vector<extension::DriftChamberDigiLocalData, std::allocator<extension::DriftChamberDigiLocalData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<extension::DriftChamberDigiLocalData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEextensioncLcLDriftChamberDigiLocalDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<extension::DriftChamberDigiLocalData>*>(nullptr))->GetClass();
      vectorlEextensioncLcLDriftChamberDigiLocalDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEextensioncLcLDriftChamberDigiLocalDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::DriftChamberDigiLocalData> : new vector<extension::DriftChamberDigiLocalData>;
   }
   static void *newArray_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::DriftChamberDigiLocalData>[nElements] : new vector<extension::DriftChamberDigiLocalData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(void *p) {
      delete (static_cast<vector<extension::DriftChamberDigiLocalData>*>(p));
   }
   static void deleteArray_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(void *p) {
      delete [] (static_cast<vector<extension::DriftChamberDigiLocalData>*>(p));
   }
   static void destruct_vectorlEextensioncLcLDriftChamberDigiLocalDatagR(void *p) {
      typedef vector<extension::DriftChamberDigiLocalData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<extension::DriftChamberDigiLocalData>

namespace ROOT {
   static TClass *vectorlEextensioncLcLDriftChamberDigiDatagR_Dictionary();
   static void vectorlEextensioncLcLDriftChamberDigiDatagR_TClassManip(TClass*);
   static void *new_vectorlEextensioncLcLDriftChamberDigiDatagR(void *p = nullptr);
   static void *newArray_vectorlEextensioncLcLDriftChamberDigiDatagR(Long_t size, void *p);
   static void delete_vectorlEextensioncLcLDriftChamberDigiDatagR(void *p);
   static void deleteArray_vectorlEextensioncLcLDriftChamberDigiDatagR(void *p);
   static void destruct_vectorlEextensioncLcLDriftChamberDigiDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<extension::DriftChamberDigiData>*)
   {
      vector<extension::DriftChamberDigiData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<extension::DriftChamberDigiData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<extension::DriftChamberDigiData>", -2, "vector", 460,
                  typeid(vector<extension::DriftChamberDigiData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEextensioncLcLDriftChamberDigiDatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<extension::DriftChamberDigiData>), alignof(vector<extension::DriftChamberDigiData>) );
      instance.SetNew(&new_vectorlEextensioncLcLDriftChamberDigiDatagR);
      instance.SetNewArray(&newArray_vectorlEextensioncLcLDriftChamberDigiDatagR);
      instance.SetDelete(&delete_vectorlEextensioncLcLDriftChamberDigiDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEextensioncLcLDriftChamberDigiDatagR);
      instance.SetDestructor(&destruct_vectorlEextensioncLcLDriftChamberDigiDatagR);
      static_assert(alignof(vector<extension::DriftChamberDigiData>::value_type) <= 4096,
          "Class with alignment strictly greater than 4096 are currently not supported in CollectionProxy. "
          "Please report this case to the developers");
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<extension::DriftChamberDigiData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<extension::DriftChamberDigiData>","std::vector<extension::DriftChamberDigiData, std::allocator<extension::DriftChamberDigiData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<extension::DriftChamberDigiData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEextensioncLcLDriftChamberDigiDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<extension::DriftChamberDigiData>*>(nullptr))->GetClass();
      vectorlEextensioncLcLDriftChamberDigiDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEextensioncLcLDriftChamberDigiDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEextensioncLcLDriftChamberDigiDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::DriftChamberDigiData> : new vector<extension::DriftChamberDigiData>;
   }
   static void *newArray_vectorlEextensioncLcLDriftChamberDigiDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<extension::DriftChamberDigiData>[nElements] : new vector<extension::DriftChamberDigiData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEextensioncLcLDriftChamberDigiDatagR(void *p) {
      delete (static_cast<vector<extension::DriftChamberDigiData>*>(p));
   }
   static void deleteArray_vectorlEextensioncLcLDriftChamberDigiDatagR(void *p) {
      delete [] (static_cast<vector<extension::DriftChamberDigiData>*>(p));
   }
   static void destruct_vectorlEextensioncLcLDriftChamberDigiDatagR(void *p) {
      typedef vector<extension::DriftChamberDigiData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<extension::DriftChamberDigiData>

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_libextensionDict() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_libextensionDict) = RecordReadRules_libextensionDict();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_libextensionDict));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_libextensionDict_Impl() {
    static const char* headers[] = {
"0",
nullptr
    };
    static const char* includePaths[] = {
"/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/WireTrackerDigi",
"/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi",
"/cvmfs/sft-nightlies.cern.ch/lcg/views/devkey-head/Mon/x86_64-el9-gcc16-opt/include",
"/cvmfs/sft-nightlies.cern.ch/lcg/latest/ROOT/6.40.04-71121/x86_64-el9-gcc16-opt/include",
"/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/build-unity-o2/WireTrackerDigi/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libextensionDict dictionary forward declarations' payload"

#pragma diagnostic push
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiData.h")))  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigi.h")))  DriftChamberDigiData;}
namespace std{template <typename _Tp> class  __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiLocalData.h")))  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiLocal.h")))  DriftChamberDigiLocalData;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiV2Data.h")))  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiV2.h")))  DriftChamberDigiV2Data;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiAssociationData.h")))  __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiAssociation.h")))  MCRecoDriftChamberDigiAssociationData;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiV2AssociationData.h")))  __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiV2Association.h")))  MCRecoDriftChamberDigiV2AssociationData;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/SenseWireHitData.h")))  __attribute__((annotate("$clingAutoload$extension/SenseWireHit.h")))  SenseWireHitData;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/SenseWireHitSimTrackerHitLinkData.h")))  __attribute__((annotate("$clingAutoload$extension/SenseWireHitSimTrackerHitLink.h")))  SenseWireHitSimTrackerHitLinkData;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/TrackData.h")))  __attribute__((annotate("$clingAutoload$extension/Track.h")))  TrackData;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigi.h")))  DriftChamberDigi;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MutableDriftChamberDigi.h")))  MutableDriftChamberDigi;}
namespace extension{class  DriftChamberDigiCollection;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiLocal.h")))  DriftChamberDigiLocal;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MutableDriftChamberDigiLocal.h")))  MutableDriftChamberDigiLocal;}
namespace extension{class  DriftChamberDigiLocalCollection;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiV2.h")))  DriftChamberDigiV2;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MutableDriftChamberDigiV2.h")))  MutableDriftChamberDigiV2;}
namespace extension{class  DriftChamberDigiV2Collection;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiAssociation.h")))  MCRecoDriftChamberDigiAssociation;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MutableMCRecoDriftChamberDigiAssociation.h")))  MutableMCRecoDriftChamberDigiAssociation;}
namespace extension{class  MCRecoDriftChamberDigiAssociationCollection;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiV2Association.h")))  MCRecoDriftChamberDigiV2Association;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MutableMCRecoDriftChamberDigiV2Association.h")))  MutableMCRecoDriftChamberDigiV2Association;}
namespace extension{class  MCRecoDriftChamberDigiV2AssociationCollection;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/SenseWireHit.h")))  SenseWireHit;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MutableSenseWireHit.h")))  MutableSenseWireHit;}
namespace extension{class  SenseWireHitCollection;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/SenseWireHitSimTrackerHitLink.h")))  SenseWireHitSimTrackerHitLink;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MutableSenseWireHitSimTrackerHitLink.h")))  MutableSenseWireHitSimTrackerHitLink;}
namespace extension{class  SenseWireHitSimTrackerHitLinkCollection;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/Track.h")))  Track;}
namespace extension{class  __attribute__((annotate("$clingAutoload$extension/MutableTrack.h")))  MutableTrack;}
namespace extension{class  TrackCollection;}
namespace extension{namespace v1{using DriftChamberDigiData __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiData.h")))  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigi.h")))  = extension::DriftChamberDigiData;}}
namespace extension{namespace v1{using DriftChamberDigiLocalData __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiLocalData.h")))  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiLocal.h")))  = extension::DriftChamberDigiLocalData;}}
namespace extension{namespace v1{using DriftChamberDigiV2Data __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiV2Data.h")))  __attribute__((annotate("$clingAutoload$extension/DriftChamberDigiV2.h")))  = extension::DriftChamberDigiV2Data;}}
namespace extension{namespace v1{using MCRecoDriftChamberDigiAssociationData __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiAssociationData.h")))  __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiAssociation.h")))  = extension::MCRecoDriftChamberDigiAssociationData;}}
namespace extension{namespace v1{using MCRecoDriftChamberDigiV2AssociationData __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiV2AssociationData.h")))  __attribute__((annotate("$clingAutoload$extension/MCRecoDriftChamberDigiV2Association.h")))  = extension::MCRecoDriftChamberDigiV2AssociationData;}}
namespace extension{namespace v1{using SenseWireHitData __attribute__((annotate("$clingAutoload$extension/SenseWireHitData.h")))  __attribute__((annotate("$clingAutoload$extension/SenseWireHit.h")))  = extension::SenseWireHitData;}}
namespace extension{namespace v1{using SenseWireHitSimTrackerHitLinkData __attribute__((annotate("$clingAutoload$extension/SenseWireHitSimTrackerHitLinkData.h")))  __attribute__((annotate("$clingAutoload$extension/SenseWireHitSimTrackerHitLink.h")))  = extension::SenseWireHitSimTrackerHitLinkData;}}
namespace extension{namespace v1{using TrackData __attribute__((annotate("$clingAutoload$extension/TrackData.h")))  __attribute__((annotate("$clingAutoload$extension/Track.h")))  = extension::TrackData;}}
#pragma diagnostic pop
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libextensionDict dictionary payload"

#ifndef PODIO_ENABLE_ROOT
  #define PODIO_ENABLE_ROOT 1
#endif
#ifndef PODIO_JSON_OUTPUT
  #define PODIO_JSON_OUTPUT 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiCollection_H
#define EXTENSION_DriftChamberDigiCollection_H

// datamodel specific includes
#include "extension/DriftChamberDigi.h"
#include "extension/MutableDriftChamberDigi.h"
#include "extension/DriftChamberDigiObj.h"
#include "extension/DriftChamberDigiCollectionData.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/detail/Pythonizations.h"
#include "podio/utilities/TypeHelpers.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <mutex>
#include <memory>
#include <cstddef>

namespace podio {
  struct RelationNames;
}

namespace extension {

/// \addtogroup extension_user_facing
/// @{

class DriftChamberDigiCollectionIterator {
public:
  using value_type = DriftChamberDigi;
  using difference_type = ptrdiff_t;
  using reference = DriftChamberDigi;
  using pointer = DriftChamberDigi*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  DriftChamberDigiCollectionIterator(size_t index, const DriftChamberDigiObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<DriftChamberDigiObj>{nullptr}), m_collection(collection) {}
  DriftChamberDigiCollectionIterator() = default;

  DriftChamberDigiCollectionIterator(const DriftChamberDigiCollectionIterator&) = default;
  DriftChamberDigiCollectionIterator(DriftChamberDigiCollectionIterator&&) = default;
  DriftChamberDigiCollectionIterator& operator=(const DriftChamberDigiCollectionIterator&) = default;
  DriftChamberDigiCollectionIterator& operator=(DriftChamberDigiCollectionIterator&&) = default;
  ~DriftChamberDigiCollectionIterator() = default;

  constexpr auto operator<=>(const DriftChamberDigiCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const DriftChamberDigiCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  DriftChamberDigiCollectionIterator& operator++();
  DriftChamberDigiCollectionIterator operator++(int);
  DriftChamberDigiCollectionIterator& operator--();
  DriftChamberDigiCollectionIterator operator--(int);
  DriftChamberDigiCollectionIterator& operator+=(difference_type n);
  DriftChamberDigiCollectionIterator operator+(difference_type n) const;
  friend DriftChamberDigiCollectionIterator operator+(difference_type n, const DriftChamberDigiCollectionIterator& it);
  DriftChamberDigiCollectionIterator& operator-=(difference_type n);
  DriftChamberDigiCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const DriftChamberDigiCollectionIterator& other) const;

private:
  size_t m_index{0};
  DriftChamberDigi m_object { podio::utils::MaybeSharedPtr<DriftChamberDigiObj>{nullptr} };
  const DriftChamberDigiObjPointerContainer* m_collection{nullptr};
};


class DriftChamberDigiMutableCollectionIterator {
public:
  using value_type = DriftChamberDigi;
  using difference_type = ptrdiff_t;
  using reference = MutableDriftChamberDigi;
  using pointer = MutableDriftChamberDigi*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  DriftChamberDigiMutableCollectionIterator(size_t index, const DriftChamberDigiObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<DriftChamberDigiObj>{nullptr}), m_collection(collection) {}
  DriftChamberDigiMutableCollectionIterator() = default;

  DriftChamberDigiMutableCollectionIterator(const DriftChamberDigiMutableCollectionIterator&) = default;
  DriftChamberDigiMutableCollectionIterator(DriftChamberDigiMutableCollectionIterator&&) = default;
  DriftChamberDigiMutableCollectionIterator& operator=(const DriftChamberDigiMutableCollectionIterator&) = default;
  DriftChamberDigiMutableCollectionIterator& operator=(DriftChamberDigiMutableCollectionIterator&&) = default;
  ~DriftChamberDigiMutableCollectionIterator() = default;

  constexpr auto operator<=>(const DriftChamberDigiMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const DriftChamberDigiMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  DriftChamberDigiMutableCollectionIterator& operator++();
  DriftChamberDigiMutableCollectionIterator operator++(int);
  DriftChamberDigiMutableCollectionIterator& operator--();
  DriftChamberDigiMutableCollectionIterator operator--(int);
  DriftChamberDigiMutableCollectionIterator& operator+=(difference_type n);
  DriftChamberDigiMutableCollectionIterator operator+(difference_type n) const;
  friend DriftChamberDigiMutableCollectionIterator operator+(difference_type n, const DriftChamberDigiMutableCollectionIterator& it);
  DriftChamberDigiMutableCollectionIterator& operator-=(difference_type n);
  DriftChamberDigiMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const DriftChamberDigiMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableDriftChamberDigi m_object { podio::utils::MaybeSharedPtr<DriftChamberDigiObj>{nullptr} };
  const DriftChamberDigiObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class DriftChamberDigiCollection : public podio::CollectionBase {
public:
  using value_type = DriftChamberDigi;
  using mutable_type = MutableDriftChamberDigi;
  using const_iterator = DriftChamberDigiCollectionIterator;
  using iterator = DriftChamberDigiMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  DriftChamberDigiCollection() = default;
  DriftChamberDigiCollection(DriftChamberDigiCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  DriftChamberDigiCollection(const DriftChamberDigiCollection& ) = delete;
  DriftChamberDigiCollection& operator=(const DriftChamberDigiCollection& ) = delete;
  DriftChamberDigiCollection(DriftChamberDigiCollection&&) = default;
  DriftChamberDigiCollection& operator=(DriftChamberDigiCollection&&) = default;

  ~DriftChamberDigiCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  DriftChamberDigiCollection(std::from_range_t, R&& range);
#endif

  /// Construct a DriftChamberDigi collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static DriftChamberDigiCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::DriftChamberDigiCollection";
  constexpr static std::string_view valueTypeName = "extension::DriftChamberDigi";
  constexpr static std::string_view dataTypeName = "extension::DriftChamberDigiData";

  /// Invalidates all iterators.
  void clear() final;

  /// Cppyy protocol to setup the pythonizations for this class. Not to be called directly.
  static void __cppjit_pythonize__(PyObject* klass, const std::string& name) {
          podio::detail::pythonizations::pythonize_subscript(klass, name);
  }
  /// Pythonizations for ROOT <= 6.40 before the switch to cppjit
  static void __cppyy_pythonize__(PyObject* klass, const std::string& name) {
        __cppjit_pythonize__(klass, name);
  }

  /// Print this collection to the passed stream
  void print(std::ostream& os=std::cout, bool flush=true) const final;

  /// Append a new object to the collection, and return this object.
  /// Does not invalidate iterators.
  MutableDriftChamberDigi create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableDriftChamberDigi create(Args&&... args);

  /// number of elements in the collection
  std::size_t size() const final;

  /// maximal number of elements in the collection
  std::size_t max_size() const final;

  /// Is the collection empty
  bool empty() const final;

  /// fully qualified type name
  const std::string_view getTypeName() const final { return typeName; }
  /// fully qualified type name of elements - with namespace
  const std::string_view getValueTypeName() const final { return valueTypeName; }
  /// fully qualified type name of stored POD elements - with namespace
  const std::string_view getDataTypeName() const final { return dataTypeName; }
  /// schema version
  podio::SchemaVersionT getSchemaVersion() const final;

  bool isSubsetCollection() const final {
    return m_isSubsetColl;
  }

  void setSubsetCollection(bool setSubset=true) final;

  /// Returns the const object of given index
  DriftChamberDigi operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableDriftChamberDigi operator[](std::size_t index);
  /// Returns the const object of given index
  DriftChamberDigi at(std::size_t index) const;
  /// Returns the object of given index
  MutableDriftChamberDigi at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableDriftChamberDigi& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const DriftChamberDigi& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (DriftChamberDigiObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
      );
    }
  }

  uint32_t getID() const final {
    return m_collectionID;
  }

  /// check if the collection has a valid ID
  bool hasID() const final {
    return getID() != static_cast<uint32_t>(podio::ObjectID::untracked) &&
        getID() != static_cast<uint32_t>(podio::ObjectID::invalid);
  }

  size_t getDatamodelRegistryIndex() const final;

  // support for the iterator protocol
  iterator begin() {
    return iterator(0, &m_storage.entries);
  }
  const_iterator begin() const {
    return const_iterator(0, &m_storage.entries);
  }
  const_iterator cbegin() const {
    return begin();
  }
  iterator end() {
    return iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator end() const {
    return const_iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator cend() const {
    return end();
  }
  // reverse iterators
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return rbegin();
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return rend();
  }


  std::vector<std::uint64_t> cellID(const size_t nElem = 0) const;
  std::vector<edm4hep::Vector3d> leftPosition(const size_t nElem = 0) const;
  std::vector<edm4hep::Vector3d> rightPosition(const size_t nElem = 0) const;
  std::vector<float> time(const size_t nElem = 0) const;
  std::vector<float> eDep(const size_t nElem = 0) const;
  std::vector<float> eDepError(const size_t nElem = 0) const;
  std::vector<std::uint32_t> clusterCount(const size_t nElem = 0) const;
  std::vector<std::int32_t> type(const size_t nElem = 0) const;
  std::vector<std::int32_t> quality(const size_t nElem = 0) const;
  std::vector<edm4hep::Vector3d> position(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class DriftChamberDigiCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable DriftChamberDigiCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const DriftChamberDigiCollection& v);

template<typename... Args>
MutableDriftChamberDigi DriftChamberDigiCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new DriftChamberDigiObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  return MutableDriftChamberDigi(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<DriftChamberDigiCollection::value_type> R>
DriftChamberDigiCollection::DriftChamberDigiCollection(std::from_range_t, R&& range)
  : DriftChamberDigiCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<DriftChamberDigiCollection::value_type> R>
DriftChamberDigiCollection DriftChamberDigiCollection::from(R&& range) {
  DriftChamberDigiCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiCollection& collection);
#endif

/// @}
} // namespace extension


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiCollection::typeName;
constexpr std::string_view extension::DriftChamberDigiCollection::valueTypeName;
constexpr std::string_view extension::DriftChamberDigiCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiCollection::typeName;
constexpr std::string_view extension::DriftChamberDigiCollection::valueTypeName;
constexpr std::string_view extension::DriftChamberDigiCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiLocalCollection_H
#define EXTENSION_DriftChamberDigiLocalCollection_H

// datamodel specific includes
#include "extension/DriftChamberDigiLocal.h"
#include "extension/MutableDriftChamberDigiLocal.h"
#include "extension/DriftChamberDigiLocalObj.h"
#include "extension/DriftChamberDigiLocalCollectionData.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/detail/Pythonizations.h"
#include "podio/utilities/TypeHelpers.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <mutex>
#include <memory>
#include <cstddef>

namespace podio {
  struct RelationNames;
}

namespace extension {

/// \addtogroup extension_user_facing
/// @{

class DriftChamberDigiLocalCollectionIterator {
public:
  using value_type = DriftChamberDigiLocal;
  using difference_type = ptrdiff_t;
  using reference = DriftChamberDigiLocal;
  using pointer = DriftChamberDigiLocal*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  DriftChamberDigiLocalCollectionIterator(size_t index, const DriftChamberDigiLocalObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>{nullptr}), m_collection(collection) {}
  DriftChamberDigiLocalCollectionIterator() = default;

  DriftChamberDigiLocalCollectionIterator(const DriftChamberDigiLocalCollectionIterator&) = default;
  DriftChamberDigiLocalCollectionIterator(DriftChamberDigiLocalCollectionIterator&&) = default;
  DriftChamberDigiLocalCollectionIterator& operator=(const DriftChamberDigiLocalCollectionIterator&) = default;
  DriftChamberDigiLocalCollectionIterator& operator=(DriftChamberDigiLocalCollectionIterator&&) = default;
  ~DriftChamberDigiLocalCollectionIterator() = default;

  constexpr auto operator<=>(const DriftChamberDigiLocalCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const DriftChamberDigiLocalCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  DriftChamberDigiLocalCollectionIterator& operator++();
  DriftChamberDigiLocalCollectionIterator operator++(int);
  DriftChamberDigiLocalCollectionIterator& operator--();
  DriftChamberDigiLocalCollectionIterator operator--(int);
  DriftChamberDigiLocalCollectionIterator& operator+=(difference_type n);
  DriftChamberDigiLocalCollectionIterator operator+(difference_type n) const;
  friend DriftChamberDigiLocalCollectionIterator operator+(difference_type n, const DriftChamberDigiLocalCollectionIterator& it);
  DriftChamberDigiLocalCollectionIterator& operator-=(difference_type n);
  DriftChamberDigiLocalCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const DriftChamberDigiLocalCollectionIterator& other) const;

private:
  size_t m_index{0};
  DriftChamberDigiLocal m_object { podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>{nullptr} };
  const DriftChamberDigiLocalObjPointerContainer* m_collection{nullptr};
};


class DriftChamberDigiLocalMutableCollectionIterator {
public:
  using value_type = DriftChamberDigiLocal;
  using difference_type = ptrdiff_t;
  using reference = MutableDriftChamberDigiLocal;
  using pointer = MutableDriftChamberDigiLocal*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  DriftChamberDigiLocalMutableCollectionIterator(size_t index, const DriftChamberDigiLocalObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>{nullptr}), m_collection(collection) {}
  DriftChamberDigiLocalMutableCollectionIterator() = default;

  DriftChamberDigiLocalMutableCollectionIterator(const DriftChamberDigiLocalMutableCollectionIterator&) = default;
  DriftChamberDigiLocalMutableCollectionIterator(DriftChamberDigiLocalMutableCollectionIterator&&) = default;
  DriftChamberDigiLocalMutableCollectionIterator& operator=(const DriftChamberDigiLocalMutableCollectionIterator&) = default;
  DriftChamberDigiLocalMutableCollectionIterator& operator=(DriftChamberDigiLocalMutableCollectionIterator&&) = default;
  ~DriftChamberDigiLocalMutableCollectionIterator() = default;

  constexpr auto operator<=>(const DriftChamberDigiLocalMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const DriftChamberDigiLocalMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  DriftChamberDigiLocalMutableCollectionIterator& operator++();
  DriftChamberDigiLocalMutableCollectionIterator operator++(int);
  DriftChamberDigiLocalMutableCollectionIterator& operator--();
  DriftChamberDigiLocalMutableCollectionIterator operator--(int);
  DriftChamberDigiLocalMutableCollectionIterator& operator+=(difference_type n);
  DriftChamberDigiLocalMutableCollectionIterator operator+(difference_type n) const;
  friend DriftChamberDigiLocalMutableCollectionIterator operator+(difference_type n, const DriftChamberDigiLocalMutableCollectionIterator& it);
  DriftChamberDigiLocalMutableCollectionIterator& operator-=(difference_type n);
  DriftChamberDigiLocalMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const DriftChamberDigiLocalMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableDriftChamberDigiLocal m_object { podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>{nullptr} };
  const DriftChamberDigiLocalObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class DriftChamberDigiLocalCollection : public podio::CollectionBase {
public:
  using value_type = DriftChamberDigiLocal;
  using mutable_type = MutableDriftChamberDigiLocal;
  using const_iterator = DriftChamberDigiLocalCollectionIterator;
  using iterator = DriftChamberDigiLocalMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  DriftChamberDigiLocalCollection() = default;
  DriftChamberDigiLocalCollection(DriftChamberDigiLocalCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  DriftChamberDigiLocalCollection(const DriftChamberDigiLocalCollection& ) = delete;
  DriftChamberDigiLocalCollection& operator=(const DriftChamberDigiLocalCollection& ) = delete;
  DriftChamberDigiLocalCollection(DriftChamberDigiLocalCollection&&) = default;
  DriftChamberDigiLocalCollection& operator=(DriftChamberDigiLocalCollection&&) = default;

  ~DriftChamberDigiLocalCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  DriftChamberDigiLocalCollection(std::from_range_t, R&& range);
#endif

  /// Construct a DriftChamberDigiLocal collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static DriftChamberDigiLocalCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::DriftChamberDigiLocalCollection";
  constexpr static std::string_view valueTypeName = "extension::DriftChamberDigiLocal";
  constexpr static std::string_view dataTypeName = "extension::DriftChamberDigiLocalData";

  /// Invalidates all iterators.
  void clear() final;

  /// Cppyy protocol to setup the pythonizations for this class. Not to be called directly.
  static void __cppjit_pythonize__(PyObject* klass, const std::string& name) {
          podio::detail::pythonizations::pythonize_subscript(klass, name);
  }
  /// Pythonizations for ROOT <= 6.40 before the switch to cppjit
  static void __cppyy_pythonize__(PyObject* klass, const std::string& name) {
        __cppjit_pythonize__(klass, name);
  }

  /// Print this collection to the passed stream
  void print(std::ostream& os=std::cout, bool flush=true) const final;

  /// Append a new object to the collection, and return this object.
  /// Does not invalidate iterators.
  MutableDriftChamberDigiLocal create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableDriftChamberDigiLocal create(Args&&... args);

  /// number of elements in the collection
  std::size_t size() const final;

  /// maximal number of elements in the collection
  std::size_t max_size() const final;

  /// Is the collection empty
  bool empty() const final;

  /// fully qualified type name
  const std::string_view getTypeName() const final { return typeName; }
  /// fully qualified type name of elements - with namespace
  const std::string_view getValueTypeName() const final { return valueTypeName; }
  /// fully qualified type name of stored POD elements - with namespace
  const std::string_view getDataTypeName() const final { return dataTypeName; }
  /// schema version
  podio::SchemaVersionT getSchemaVersion() const final;

  bool isSubsetCollection() const final {
    return m_isSubsetColl;
  }

  void setSubsetCollection(bool setSubset=true) final;

  /// Returns the const object of given index
  DriftChamberDigiLocal operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableDriftChamberDigiLocal operator[](std::size_t index);
  /// Returns the const object of given index
  DriftChamberDigiLocal at(std::size_t index) const;
  /// Returns the object of given index
  MutableDriftChamberDigiLocal at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableDriftChamberDigiLocal& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const DriftChamberDigiLocal& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (DriftChamberDigiLocalObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
      );
    }
  }

  uint32_t getID() const final {
    return m_collectionID;
  }

  /// check if the collection has a valid ID
  bool hasID() const final {
    return getID() != static_cast<uint32_t>(podio::ObjectID::untracked) &&
        getID() != static_cast<uint32_t>(podio::ObjectID::invalid);
  }

  size_t getDatamodelRegistryIndex() const final;

  // support for the iterator protocol
  iterator begin() {
    return iterator(0, &m_storage.entries);
  }
  const_iterator begin() const {
    return const_iterator(0, &m_storage.entries);
  }
  const_iterator cbegin() const {
    return begin();
  }
  iterator end() {
    return iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator end() const {
    return const_iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator cend() const {
    return end();
  }
  // reverse iterators
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return rbegin();
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return rend();
  }


  std::vector<std::uint64_t> cellID(const size_t nElem = 0) const;
  std::vector<float> distanceToWire(const size_t nElem = 0) const;
  std::vector<float> zPositionAlongWire(const size_t nElem = 0) const;
  std::vector<float> time(const size_t nElem = 0) const;
  std::vector<float> eDep(const size_t nElem = 0) const;
  std::vector<float> eDepError(const size_t nElem = 0) const;
  std::vector<std::uint32_t> clusterCount(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class DriftChamberDigiLocalCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable DriftChamberDigiLocalCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const DriftChamberDigiLocalCollection& v);

template<typename... Args>
MutableDriftChamberDigiLocal DriftChamberDigiLocalCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new DriftChamberDigiLocalObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  return MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<DriftChamberDigiLocalCollection::value_type> R>
DriftChamberDigiLocalCollection::DriftChamberDigiLocalCollection(std::from_range_t, R&& range)
  : DriftChamberDigiLocalCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<DriftChamberDigiLocalCollection::value_type> R>
DriftChamberDigiLocalCollection DriftChamberDigiLocalCollection::from(R&& range) {
  DriftChamberDigiLocalCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiLocalCollection& collection);
#endif

/// @}
} // namespace extension


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiLocalCollection::typeName;
constexpr std::string_view extension::DriftChamberDigiLocalCollection::valueTypeName;
constexpr std::string_view extension::DriftChamberDigiLocalCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiLocalCollection::typeName;
constexpr std::string_view extension::DriftChamberDigiLocalCollection::valueTypeName;
constexpr std::string_view extension::DriftChamberDigiLocalCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2Collection_H
#define EXTENSION_DriftChamberDigiV2Collection_H

// datamodel specific includes
#include "extension/DriftChamberDigiV2.h"
#include "extension/MutableDriftChamberDigiV2.h"
#include "extension/DriftChamberDigiV2Obj.h"
#include "extension/DriftChamberDigiV2CollectionData.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/detail/Pythonizations.h"
#include "podio/utilities/TypeHelpers.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <mutex>
#include <memory>
#include <cstddef>

namespace podio {
  struct RelationNames;
}

namespace extension {

/// \addtogroup extension_user_facing
/// @{

class DriftChamberDigiV2CollectionIterator {
public:
  using value_type = DriftChamberDigiV2;
  using difference_type = ptrdiff_t;
  using reference = DriftChamberDigiV2;
  using pointer = DriftChamberDigiV2*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  DriftChamberDigiV2CollectionIterator(size_t index, const DriftChamberDigiV2ObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr}), m_collection(collection) {}
  DriftChamberDigiV2CollectionIterator() = default;

  DriftChamberDigiV2CollectionIterator(const DriftChamberDigiV2CollectionIterator&) = default;
  DriftChamberDigiV2CollectionIterator(DriftChamberDigiV2CollectionIterator&&) = default;
  DriftChamberDigiV2CollectionIterator& operator=(const DriftChamberDigiV2CollectionIterator&) = default;
  DriftChamberDigiV2CollectionIterator& operator=(DriftChamberDigiV2CollectionIterator&&) = default;
  ~DriftChamberDigiV2CollectionIterator() = default;

  constexpr auto operator<=>(const DriftChamberDigiV2CollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const DriftChamberDigiV2CollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  DriftChamberDigiV2CollectionIterator& operator++();
  DriftChamberDigiV2CollectionIterator operator++(int);
  DriftChamberDigiV2CollectionIterator& operator--();
  DriftChamberDigiV2CollectionIterator operator--(int);
  DriftChamberDigiV2CollectionIterator& operator+=(difference_type n);
  DriftChamberDigiV2CollectionIterator operator+(difference_type n) const;
  friend DriftChamberDigiV2CollectionIterator operator+(difference_type n, const DriftChamberDigiV2CollectionIterator& it);
  DriftChamberDigiV2CollectionIterator& operator-=(difference_type n);
  DriftChamberDigiV2CollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const DriftChamberDigiV2CollectionIterator& other) const;

private:
  size_t m_index{0};
  DriftChamberDigiV2 m_object { podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr} };
  const DriftChamberDigiV2ObjPointerContainer* m_collection{nullptr};
};


class DriftChamberDigiV2MutableCollectionIterator {
public:
  using value_type = DriftChamberDigiV2;
  using difference_type = ptrdiff_t;
  using reference = MutableDriftChamberDigiV2;
  using pointer = MutableDriftChamberDigiV2*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  DriftChamberDigiV2MutableCollectionIterator(size_t index, const DriftChamberDigiV2ObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr}), m_collection(collection) {}
  DriftChamberDigiV2MutableCollectionIterator() = default;

  DriftChamberDigiV2MutableCollectionIterator(const DriftChamberDigiV2MutableCollectionIterator&) = default;
  DriftChamberDigiV2MutableCollectionIterator(DriftChamberDigiV2MutableCollectionIterator&&) = default;
  DriftChamberDigiV2MutableCollectionIterator& operator=(const DriftChamberDigiV2MutableCollectionIterator&) = default;
  DriftChamberDigiV2MutableCollectionIterator& operator=(DriftChamberDigiV2MutableCollectionIterator&&) = default;
  ~DriftChamberDigiV2MutableCollectionIterator() = default;

  constexpr auto operator<=>(const DriftChamberDigiV2MutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const DriftChamberDigiV2MutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  DriftChamberDigiV2MutableCollectionIterator& operator++();
  DriftChamberDigiV2MutableCollectionIterator operator++(int);
  DriftChamberDigiV2MutableCollectionIterator& operator--();
  DriftChamberDigiV2MutableCollectionIterator operator--(int);
  DriftChamberDigiV2MutableCollectionIterator& operator+=(difference_type n);
  DriftChamberDigiV2MutableCollectionIterator operator+(difference_type n) const;
  friend DriftChamberDigiV2MutableCollectionIterator operator+(difference_type n, const DriftChamberDigiV2MutableCollectionIterator& it);
  DriftChamberDigiV2MutableCollectionIterator& operator-=(difference_type n);
  DriftChamberDigiV2MutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const DriftChamberDigiV2MutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableDriftChamberDigiV2 m_object { podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr} };
  const DriftChamberDigiV2ObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class DriftChamberDigiV2Collection : public podio::CollectionBase {
public:
  using value_type = DriftChamberDigiV2;
  using mutable_type = MutableDriftChamberDigiV2;
  using const_iterator = DriftChamberDigiV2CollectionIterator;
  using iterator = DriftChamberDigiV2MutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  DriftChamberDigiV2Collection() = default;
  DriftChamberDigiV2Collection(DriftChamberDigiV2CollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  DriftChamberDigiV2Collection(const DriftChamberDigiV2Collection& ) = delete;
  DriftChamberDigiV2Collection& operator=(const DriftChamberDigiV2Collection& ) = delete;
  DriftChamberDigiV2Collection(DriftChamberDigiV2Collection&&) = default;
  DriftChamberDigiV2Collection& operator=(DriftChamberDigiV2Collection&&) = default;

  ~DriftChamberDigiV2Collection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  DriftChamberDigiV2Collection(std::from_range_t, R&& range);
#endif

  /// Construct a DriftChamberDigiV2 collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static DriftChamberDigiV2Collection from(R&& range);

  constexpr static std::string_view typeName = "extension::DriftChamberDigiV2Collection";
  constexpr static std::string_view valueTypeName = "extension::DriftChamberDigiV2";
  constexpr static std::string_view dataTypeName = "extension::DriftChamberDigiV2Data";

  /// Invalidates all iterators.
  void clear() final;

  /// Cppyy protocol to setup the pythonizations for this class. Not to be called directly.
  static void __cppjit_pythonize__(PyObject* klass, const std::string& name) {
          podio::detail::pythonizations::pythonize_subscript(klass, name);
  }
  /// Pythonizations for ROOT <= 6.40 before the switch to cppjit
  static void __cppyy_pythonize__(PyObject* klass, const std::string& name) {
        __cppjit_pythonize__(klass, name);
  }

  /// Print this collection to the passed stream
  void print(std::ostream& os=std::cout, bool flush=true) const final;

  /// Append a new object to the collection, and return this object.
  /// Does not invalidate iterators.
  MutableDriftChamberDigiV2 create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableDriftChamberDigiV2 create(Args&&... args);

  /// number of elements in the collection
  std::size_t size() const final;

  /// maximal number of elements in the collection
  std::size_t max_size() const final;

  /// Is the collection empty
  bool empty() const final;

  /// fully qualified type name
  const std::string_view getTypeName() const final { return typeName; }
  /// fully qualified type name of elements - with namespace
  const std::string_view getValueTypeName() const final { return valueTypeName; }
  /// fully qualified type name of stored POD elements - with namespace
  const std::string_view getDataTypeName() const final { return dataTypeName; }
  /// schema version
  podio::SchemaVersionT getSchemaVersion() const final;

  bool isSubsetCollection() const final {
    return m_isSubsetColl;
  }

  void setSubsetCollection(bool setSubset=true) final;

  /// Returns the const object of given index
  DriftChamberDigiV2 operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableDriftChamberDigiV2 operator[](std::size_t index);
  /// Returns the const object of given index
  DriftChamberDigiV2 at(std::size_t index) const;
  /// Returns the object of given index
  MutableDriftChamberDigiV2 at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableDriftChamberDigiV2& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const DriftChamberDigiV2& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (DriftChamberDigiV2Obj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
      );
    }
  }

  uint32_t getID() const final {
    return m_collectionID;
  }

  /// check if the collection has a valid ID
  bool hasID() const final {
    return getID() != static_cast<uint32_t>(podio::ObjectID::untracked) &&
        getID() != static_cast<uint32_t>(podio::ObjectID::invalid);
  }

  size_t getDatamodelRegistryIndex() const final;

  // support for the iterator protocol
  iterator begin() {
    return iterator(0, &m_storage.entries);
  }
  const_iterator begin() const {
    return const_iterator(0, &m_storage.entries);
  }
  const_iterator cbegin() const {
    return begin();
  }
  iterator end() {
    return iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator end() const {
    return const_iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator cend() const {
    return end();
  }
  // reverse iterators
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return rbegin();
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return rend();
  }


  std::vector<std::uint64_t> cellID(const size_t nElem = 0) const;
  std::vector<std::int32_t> type(const size_t nElem = 0) const;
  std::vector<std::int32_t> quality(const size_t nElem = 0) const;
  std::vector<float> time(const size_t nElem = 0) const;
  std::vector<float> eDep(const size_t nElem = 0) const;
  std::vector<float> eDepError(const size_t nElem = 0) const;
  std::vector<edm4hep::Vector3d> position(const size_t nElem = 0) const;
  std::vector<edm4hep::Vector3d> directionSW(const size_t nElem = 0) const;
  std::vector<float> distanceToWire(const size_t nElem = 0) const;
  std::vector<std::uint32_t> nCluster(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class DriftChamberDigiV2CollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable DriftChamberDigiV2CollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const DriftChamberDigiV2Collection& v);

template<typename... Args>
MutableDriftChamberDigiV2 DriftChamberDigiV2Collection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new DriftChamberDigiV2Obj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  // Need to initialize the relation vectors manually for the {ObjectID, DriftChamberDigiV2Data} constructor
  obj->m_nElectrons = new std::vector<std::uint16_t>();
  m_storage.createRelations(obj);
  return MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<DriftChamberDigiV2Collection::value_type> R>
DriftChamberDigiV2Collection::DriftChamberDigiV2Collection(std::from_range_t, R&& range)
  : DriftChamberDigiV2Collection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<DriftChamberDigiV2Collection::value_type> R>
DriftChamberDigiV2Collection DriftChamberDigiV2Collection::from(R&& range) {
  DriftChamberDigiV2Collection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiV2Collection& collection);
#endif

/// @}
} // namespace extension


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiV2Collection::typeName;
constexpr std::string_view extension::DriftChamberDigiV2Collection::valueTypeName;
constexpr std::string_view extension::DriftChamberDigiV2Collection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiV2Collection::typeName;
constexpr std::string_view extension::DriftChamberDigiV2Collection::valueTypeName;
constexpr std::string_view extension::DriftChamberDigiV2Collection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociationCollection_H
#define EXTENSION_MCRecoDriftChamberDigiAssociationCollection_H

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiAssociation.h"
#include "extension/MutableMCRecoDriftChamberDigiAssociation.h"
#include "extension/MCRecoDriftChamberDigiAssociationObj.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollectionData.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/detail/Pythonizations.h"
#include "podio/utilities/TypeHelpers.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <mutex>
#include <memory>
#include <cstddef>

namespace podio {
  struct RelationNames;
}

namespace extension {

/// \addtogroup extension_user_facing
/// @{

class MCRecoDriftChamberDigiAssociationCollectionIterator {
public:
  using value_type = MCRecoDriftChamberDigiAssociation;
  using difference_type = ptrdiff_t;
  using reference = MCRecoDriftChamberDigiAssociation;
  using pointer = MCRecoDriftChamberDigiAssociation*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  MCRecoDriftChamberDigiAssociationCollectionIterator(size_t index, const MCRecoDriftChamberDigiAssociationObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr}), m_collection(collection) {}
  MCRecoDriftChamberDigiAssociationCollectionIterator() = default;

  MCRecoDriftChamberDigiAssociationCollectionIterator(const MCRecoDriftChamberDigiAssociationCollectionIterator&) = default;
  MCRecoDriftChamberDigiAssociationCollectionIterator(MCRecoDriftChamberDigiAssociationCollectionIterator&&) = default;
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator=(const MCRecoDriftChamberDigiAssociationCollectionIterator&) = default;
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator=(MCRecoDriftChamberDigiAssociationCollectionIterator&&) = default;
  ~MCRecoDriftChamberDigiAssociationCollectionIterator() = default;

  constexpr auto operator<=>(const MCRecoDriftChamberDigiAssociationCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const MCRecoDriftChamberDigiAssociationCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator++();
  MCRecoDriftChamberDigiAssociationCollectionIterator operator++(int);
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator--();
  MCRecoDriftChamberDigiAssociationCollectionIterator operator--(int);
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator+=(difference_type n);
  MCRecoDriftChamberDigiAssociationCollectionIterator operator+(difference_type n) const;
  friend MCRecoDriftChamberDigiAssociationCollectionIterator operator+(difference_type n, const MCRecoDriftChamberDigiAssociationCollectionIterator& it);
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator-=(difference_type n);
  MCRecoDriftChamberDigiAssociationCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const MCRecoDriftChamberDigiAssociationCollectionIterator& other) const;

private:
  size_t m_index{0};
  MCRecoDriftChamberDigiAssociation m_object { podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr} };
  const MCRecoDriftChamberDigiAssociationObjPointerContainer* m_collection{nullptr};
};


class MCRecoDriftChamberDigiAssociationMutableCollectionIterator {
public:
  using value_type = MCRecoDriftChamberDigiAssociation;
  using difference_type = ptrdiff_t;
  using reference = MutableMCRecoDriftChamberDigiAssociation;
  using pointer = MutableMCRecoDriftChamberDigiAssociation*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  MCRecoDriftChamberDigiAssociationMutableCollectionIterator(size_t index, const MCRecoDriftChamberDigiAssociationObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr}), m_collection(collection) {}
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator() = default;

  MCRecoDriftChamberDigiAssociationMutableCollectionIterator(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator&) = default;
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator(MCRecoDriftChamberDigiAssociationMutableCollectionIterator&&) = default;
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator=(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator&) = default;
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator=(MCRecoDriftChamberDigiAssociationMutableCollectionIterator&&) = default;
  ~MCRecoDriftChamberDigiAssociationMutableCollectionIterator() = default;

  constexpr auto operator<=>(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator++();
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator++(int);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator--();
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator--(int);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator+=(difference_type n);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator+(difference_type n) const;
  friend MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator+(difference_type n, const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& it);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator-=(difference_type n);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableMCRecoDriftChamberDigiAssociation m_object { podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr} };
  const MCRecoDriftChamberDigiAssociationObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class MCRecoDriftChamberDigiAssociationCollection : public podio::CollectionBase {
public:
  using value_type = MCRecoDriftChamberDigiAssociation;
  using mutable_type = MutableMCRecoDriftChamberDigiAssociation;
  using const_iterator = MCRecoDriftChamberDigiAssociationCollectionIterator;
  using iterator = MCRecoDriftChamberDigiAssociationMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  MCRecoDriftChamberDigiAssociationCollection() = default;
  MCRecoDriftChamberDigiAssociationCollection(MCRecoDriftChamberDigiAssociationCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  MCRecoDriftChamberDigiAssociationCollection(const MCRecoDriftChamberDigiAssociationCollection& ) = delete;
  MCRecoDriftChamberDigiAssociationCollection& operator=(const MCRecoDriftChamberDigiAssociationCollection& ) = delete;
  MCRecoDriftChamberDigiAssociationCollection(MCRecoDriftChamberDigiAssociationCollection&&) = default;
  MCRecoDriftChamberDigiAssociationCollection& operator=(MCRecoDriftChamberDigiAssociationCollection&&) = default;

  ~MCRecoDriftChamberDigiAssociationCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  MCRecoDriftChamberDigiAssociationCollection(std::from_range_t, R&& range);
#endif

  /// Construct a MCRecoDriftChamberDigiAssociation collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static MCRecoDriftChamberDigiAssociationCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::MCRecoDriftChamberDigiAssociationCollection";
  constexpr static std::string_view valueTypeName = "extension::MCRecoDriftChamberDigiAssociation";
  constexpr static std::string_view dataTypeName = "extension::MCRecoDriftChamberDigiAssociationData";

  /// Invalidates all iterators.
  void clear() final;

  /// Cppyy protocol to setup the pythonizations for this class. Not to be called directly.
  static void __cppjit_pythonize__(PyObject* klass, const std::string& name) {
          podio::detail::pythonizations::pythonize_subscript(klass, name);
  }
  /// Pythonizations for ROOT <= 6.40 before the switch to cppjit
  static void __cppyy_pythonize__(PyObject* klass, const std::string& name) {
        __cppjit_pythonize__(klass, name);
  }

  /// Print this collection to the passed stream
  void print(std::ostream& os=std::cout, bool flush=true) const final;

  /// Append a new object to the collection, and return this object.
  /// Does not invalidate iterators.
  MutableMCRecoDriftChamberDigiAssociation create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableMCRecoDriftChamberDigiAssociation create(Args&&... args);

  /// number of elements in the collection
  std::size_t size() const final;

  /// maximal number of elements in the collection
  std::size_t max_size() const final;

  /// Is the collection empty
  bool empty() const final;

  /// fully qualified type name
  const std::string_view getTypeName() const final { return typeName; }
  /// fully qualified type name of elements - with namespace
  const std::string_view getValueTypeName() const final { return valueTypeName; }
  /// fully qualified type name of stored POD elements - with namespace
  const std::string_view getDataTypeName() const final { return dataTypeName; }
  /// schema version
  podio::SchemaVersionT getSchemaVersion() const final;

  bool isSubsetCollection() const final {
    return m_isSubsetColl;
  }

  void setSubsetCollection(bool setSubset=true) final;

  /// Returns the const object of given index
  MCRecoDriftChamberDigiAssociation operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableMCRecoDriftChamberDigiAssociation operator[](std::size_t index);
  /// Returns the const object of given index
  MCRecoDriftChamberDigiAssociation at(std::size_t index) const;
  /// Returns the object of given index
  MutableMCRecoDriftChamberDigiAssociation at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableMCRecoDriftChamberDigiAssociation& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const MCRecoDriftChamberDigiAssociation& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (MCRecoDriftChamberDigiAssociationObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
      );
    }
  }

  uint32_t getID() const final {
    return m_collectionID;
  }

  /// check if the collection has a valid ID
  bool hasID() const final {
    return getID() != static_cast<uint32_t>(podio::ObjectID::untracked) &&
        getID() != static_cast<uint32_t>(podio::ObjectID::invalid);
  }

  size_t getDatamodelRegistryIndex() const final;

  // support for the iterator protocol
  iterator begin() {
    return iterator(0, &m_storage.entries);
  }
  const_iterator begin() const {
    return const_iterator(0, &m_storage.entries);
  }
  const_iterator cbegin() const {
    return begin();
  }
  iterator end() {
    return iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator end() const {
    return const_iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator cend() const {
    return end();
  }
  // reverse iterators
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return rbegin();
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return rend();
  }


  std::vector<float> weight(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class MCRecoDriftChamberDigiAssociationCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable MCRecoDriftChamberDigiAssociationCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiAssociationCollection& v);

template<typename... Args>
MutableMCRecoDriftChamberDigiAssociation MCRecoDriftChamberDigiAssociationCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new MCRecoDriftChamberDigiAssociationObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  return MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<MCRecoDriftChamberDigiAssociationCollection::value_type> R>
MCRecoDriftChamberDigiAssociationCollection::MCRecoDriftChamberDigiAssociationCollection(std::from_range_t, R&& range)
  : MCRecoDriftChamberDigiAssociationCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<MCRecoDriftChamberDigiAssociationCollection::value_type> R>
MCRecoDriftChamberDigiAssociationCollection MCRecoDriftChamberDigiAssociationCollection::from(R&& range) {
  MCRecoDriftChamberDigiAssociationCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiAssociationCollection& collection);
#endif

/// @}
} // namespace extension


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::typeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::valueTypeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::typeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::valueTypeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiV2AssociationCollection_H
#define EXTENSION_MCRecoDriftChamberDigiV2AssociationCollection_H

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiV2Association.h"
#include "extension/MutableMCRecoDriftChamberDigiV2Association.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationCollectionData.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/detail/Pythonizations.h"
#include "podio/utilities/TypeHelpers.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <mutex>
#include <memory>
#include <cstddef>

namespace podio {
  struct RelationNames;
}

namespace extension {

/// \addtogroup extension_user_facing
/// @{

class MCRecoDriftChamberDigiV2AssociationCollectionIterator {
public:
  using value_type = MCRecoDriftChamberDigiV2Association;
  using difference_type = ptrdiff_t;
  using reference = MCRecoDriftChamberDigiV2Association;
  using pointer = MCRecoDriftChamberDigiV2Association*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  MCRecoDriftChamberDigiV2AssociationCollectionIterator(size_t index, const MCRecoDriftChamberDigiV2AssociationObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr}), m_collection(collection) {}
  MCRecoDriftChamberDigiV2AssociationCollectionIterator() = default;

  MCRecoDriftChamberDigiV2AssociationCollectionIterator(const MCRecoDriftChamberDigiV2AssociationCollectionIterator&) = default;
  MCRecoDriftChamberDigiV2AssociationCollectionIterator(MCRecoDriftChamberDigiV2AssociationCollectionIterator&&) = default;
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator=(const MCRecoDriftChamberDigiV2AssociationCollectionIterator&) = default;
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator=(MCRecoDriftChamberDigiV2AssociationCollectionIterator&&) = default;
  ~MCRecoDriftChamberDigiV2AssociationCollectionIterator() = default;

  constexpr auto operator<=>(const MCRecoDriftChamberDigiV2AssociationCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const MCRecoDriftChamberDigiV2AssociationCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator++();
  MCRecoDriftChamberDigiV2AssociationCollectionIterator operator++(int);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator--();
  MCRecoDriftChamberDigiV2AssociationCollectionIterator operator--(int);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator+=(difference_type n);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator operator+(difference_type n) const;
  friend MCRecoDriftChamberDigiV2AssociationCollectionIterator operator+(difference_type n, const MCRecoDriftChamberDigiV2AssociationCollectionIterator& it);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator-=(difference_type n);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const MCRecoDriftChamberDigiV2AssociationCollectionIterator& other) const;

private:
  size_t m_index{0};
  MCRecoDriftChamberDigiV2Association m_object { podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr} };
  const MCRecoDriftChamberDigiV2AssociationObjPointerContainer* m_collection{nullptr};
};


class MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator {
public:
  using value_type = MCRecoDriftChamberDigiV2Association;
  using difference_type = ptrdiff_t;
  using reference = MutableMCRecoDriftChamberDigiV2Association;
  using pointer = MutableMCRecoDriftChamberDigiV2Association*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator(size_t index, const MCRecoDriftChamberDigiV2AssociationObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr}), m_collection(collection) {}
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator() = default;

  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator&) = default;
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator(MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator&&) = default;
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator=(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator&) = default;
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator=(MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator&&) = default;
  ~MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator() = default;

  constexpr auto operator<=>(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator++();
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator++(int);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator--();
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator--(int);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator+=(difference_type n);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator+(difference_type n) const;
  friend MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator+(difference_type n, const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& it);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator-=(difference_type n);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableMCRecoDriftChamberDigiV2Association m_object { podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr} };
  const MCRecoDriftChamberDigiV2AssociationObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class MCRecoDriftChamberDigiV2AssociationCollection : public podio::CollectionBase {
public:
  using value_type = MCRecoDriftChamberDigiV2Association;
  using mutable_type = MutableMCRecoDriftChamberDigiV2Association;
  using const_iterator = MCRecoDriftChamberDigiV2AssociationCollectionIterator;
  using iterator = MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  MCRecoDriftChamberDigiV2AssociationCollection() = default;
  MCRecoDriftChamberDigiV2AssociationCollection(MCRecoDriftChamberDigiV2AssociationCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  MCRecoDriftChamberDigiV2AssociationCollection(const MCRecoDriftChamberDigiV2AssociationCollection& ) = delete;
  MCRecoDriftChamberDigiV2AssociationCollection& operator=(const MCRecoDriftChamberDigiV2AssociationCollection& ) = delete;
  MCRecoDriftChamberDigiV2AssociationCollection(MCRecoDriftChamberDigiV2AssociationCollection&&) = default;
  MCRecoDriftChamberDigiV2AssociationCollection& operator=(MCRecoDriftChamberDigiV2AssociationCollection&&) = default;

  ~MCRecoDriftChamberDigiV2AssociationCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  MCRecoDriftChamberDigiV2AssociationCollection(std::from_range_t, R&& range);
#endif

  /// Construct a MCRecoDriftChamberDigiV2Association collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static MCRecoDriftChamberDigiV2AssociationCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::MCRecoDriftChamberDigiV2AssociationCollection";
  constexpr static std::string_view valueTypeName = "extension::MCRecoDriftChamberDigiV2Association";
  constexpr static std::string_view dataTypeName = "extension::MCRecoDriftChamberDigiV2AssociationData";

  /// Invalidates all iterators.
  void clear() final;

  /// Cppyy protocol to setup the pythonizations for this class. Not to be called directly.
  static void __cppjit_pythonize__(PyObject* klass, const std::string& name) {
          podio::detail::pythonizations::pythonize_subscript(klass, name);
  }
  /// Pythonizations for ROOT <= 6.40 before the switch to cppjit
  static void __cppyy_pythonize__(PyObject* klass, const std::string& name) {
        __cppjit_pythonize__(klass, name);
  }

  /// Print this collection to the passed stream
  void print(std::ostream& os=std::cout, bool flush=true) const final;

  /// Append a new object to the collection, and return this object.
  /// Does not invalidate iterators.
  MutableMCRecoDriftChamberDigiV2Association create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableMCRecoDriftChamberDigiV2Association create(Args&&... args);

  /// number of elements in the collection
  std::size_t size() const final;

  /// maximal number of elements in the collection
  std::size_t max_size() const final;

  /// Is the collection empty
  bool empty() const final;

  /// fully qualified type name
  const std::string_view getTypeName() const final { return typeName; }
  /// fully qualified type name of elements - with namespace
  const std::string_view getValueTypeName() const final { return valueTypeName; }
  /// fully qualified type name of stored POD elements - with namespace
  const std::string_view getDataTypeName() const final { return dataTypeName; }
  /// schema version
  podio::SchemaVersionT getSchemaVersion() const final;

  bool isSubsetCollection() const final {
    return m_isSubsetColl;
  }

  void setSubsetCollection(bool setSubset=true) final;

  /// Returns the const object of given index
  MCRecoDriftChamberDigiV2Association operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableMCRecoDriftChamberDigiV2Association operator[](std::size_t index);
  /// Returns the const object of given index
  MCRecoDriftChamberDigiV2Association at(std::size_t index) const;
  /// Returns the object of given index
  MutableMCRecoDriftChamberDigiV2Association at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableMCRecoDriftChamberDigiV2Association& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const MCRecoDriftChamberDigiV2Association& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (MCRecoDriftChamberDigiV2AssociationObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
      );
    }
  }

  uint32_t getID() const final {
    return m_collectionID;
  }

  /// check if the collection has a valid ID
  bool hasID() const final {
    return getID() != static_cast<uint32_t>(podio::ObjectID::untracked) &&
        getID() != static_cast<uint32_t>(podio::ObjectID::invalid);
  }

  size_t getDatamodelRegistryIndex() const final;

  // support for the iterator protocol
  iterator begin() {
    return iterator(0, &m_storage.entries);
  }
  const_iterator begin() const {
    return const_iterator(0, &m_storage.entries);
  }
  const_iterator cbegin() const {
    return begin();
  }
  iterator end() {
    return iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator end() const {
    return const_iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator cend() const {
    return end();
  }
  // reverse iterators
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return rbegin();
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return rend();
  }


  std::vector<float> weight(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class MCRecoDriftChamberDigiV2AssociationCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable MCRecoDriftChamberDigiV2AssociationCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiV2AssociationCollection& v);

template<typename... Args>
MutableMCRecoDriftChamberDigiV2Association MCRecoDriftChamberDigiV2AssociationCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new MCRecoDriftChamberDigiV2AssociationObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  return MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<MCRecoDriftChamberDigiV2AssociationCollection::value_type> R>
MCRecoDriftChamberDigiV2AssociationCollection::MCRecoDriftChamberDigiV2AssociationCollection(std::from_range_t, R&& range)
  : MCRecoDriftChamberDigiV2AssociationCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<MCRecoDriftChamberDigiV2AssociationCollection::value_type> R>
MCRecoDriftChamberDigiV2AssociationCollection MCRecoDriftChamberDigiV2AssociationCollection::from(R&& range) {
  MCRecoDriftChamberDigiV2AssociationCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiV2AssociationCollection& collection);
#endif

/// @}
} // namespace extension


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::typeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::valueTypeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::typeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::valueTypeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitCollection_H
#define EXTENSION_SenseWireHitCollection_H

// datamodel specific includes
#include "extension/SenseWireHit.h"
#include "extension/MutableSenseWireHit.h"
#include "extension/SenseWireHitObj.h"
#include "extension/SenseWireHitCollectionData.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/detail/Pythonizations.h"
#include "podio/utilities/TypeHelpers.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <mutex>
#include <memory>
#include <cstddef>

namespace podio {
  struct RelationNames;
}

namespace extension {

/// \addtogroup extension_user_facing
/// @{

class SenseWireHitCollectionIterator {
public:
  using value_type = SenseWireHit;
  using difference_type = ptrdiff_t;
  using reference = SenseWireHit;
  using pointer = SenseWireHit*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  SenseWireHitCollectionIterator(size_t index, const SenseWireHitObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr}), m_collection(collection) {}
  SenseWireHitCollectionIterator() = default;

  SenseWireHitCollectionIterator(const SenseWireHitCollectionIterator&) = default;
  SenseWireHitCollectionIterator(SenseWireHitCollectionIterator&&) = default;
  SenseWireHitCollectionIterator& operator=(const SenseWireHitCollectionIterator&) = default;
  SenseWireHitCollectionIterator& operator=(SenseWireHitCollectionIterator&&) = default;
  ~SenseWireHitCollectionIterator() = default;

  constexpr auto operator<=>(const SenseWireHitCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const SenseWireHitCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  SenseWireHitCollectionIterator& operator++();
  SenseWireHitCollectionIterator operator++(int);
  SenseWireHitCollectionIterator& operator--();
  SenseWireHitCollectionIterator operator--(int);
  SenseWireHitCollectionIterator& operator+=(difference_type n);
  SenseWireHitCollectionIterator operator+(difference_type n) const;
  friend SenseWireHitCollectionIterator operator+(difference_type n, const SenseWireHitCollectionIterator& it);
  SenseWireHitCollectionIterator& operator-=(difference_type n);
  SenseWireHitCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const SenseWireHitCollectionIterator& other) const;

private:
  size_t m_index{0};
  SenseWireHit m_object { podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr} };
  const SenseWireHitObjPointerContainer* m_collection{nullptr};
};


class SenseWireHitMutableCollectionIterator {
public:
  using value_type = SenseWireHit;
  using difference_type = ptrdiff_t;
  using reference = MutableSenseWireHit;
  using pointer = MutableSenseWireHit*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  SenseWireHitMutableCollectionIterator(size_t index, const SenseWireHitObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr}), m_collection(collection) {}
  SenseWireHitMutableCollectionIterator() = default;

  SenseWireHitMutableCollectionIterator(const SenseWireHitMutableCollectionIterator&) = default;
  SenseWireHitMutableCollectionIterator(SenseWireHitMutableCollectionIterator&&) = default;
  SenseWireHitMutableCollectionIterator& operator=(const SenseWireHitMutableCollectionIterator&) = default;
  SenseWireHitMutableCollectionIterator& operator=(SenseWireHitMutableCollectionIterator&&) = default;
  ~SenseWireHitMutableCollectionIterator() = default;

  constexpr auto operator<=>(const SenseWireHitMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const SenseWireHitMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  SenseWireHitMutableCollectionIterator& operator++();
  SenseWireHitMutableCollectionIterator operator++(int);
  SenseWireHitMutableCollectionIterator& operator--();
  SenseWireHitMutableCollectionIterator operator--(int);
  SenseWireHitMutableCollectionIterator& operator+=(difference_type n);
  SenseWireHitMutableCollectionIterator operator+(difference_type n) const;
  friend SenseWireHitMutableCollectionIterator operator+(difference_type n, const SenseWireHitMutableCollectionIterator& it);
  SenseWireHitMutableCollectionIterator& operator-=(difference_type n);
  SenseWireHitMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const SenseWireHitMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableSenseWireHit m_object { podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr} };
  const SenseWireHitObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class SenseWireHitCollection : public podio::CollectionBase {
public:
  using value_type = SenseWireHit;
  using mutable_type = MutableSenseWireHit;
  using const_iterator = SenseWireHitCollectionIterator;
  using iterator = SenseWireHitMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  SenseWireHitCollection() = default;
  SenseWireHitCollection(SenseWireHitCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  SenseWireHitCollection(const SenseWireHitCollection& ) = delete;
  SenseWireHitCollection& operator=(const SenseWireHitCollection& ) = delete;
  SenseWireHitCollection(SenseWireHitCollection&&) = default;
  SenseWireHitCollection& operator=(SenseWireHitCollection&&) = default;

  ~SenseWireHitCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  SenseWireHitCollection(std::from_range_t, R&& range);
#endif

  /// Construct a SenseWireHit collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static SenseWireHitCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::SenseWireHitCollection";
  constexpr static std::string_view valueTypeName = "extension::SenseWireHit";
  constexpr static std::string_view dataTypeName = "extension::SenseWireHitData";

  /// Invalidates all iterators.
  void clear() final;

  /// Cppyy protocol to setup the pythonizations for this class. Not to be called directly.
  static void __cppjit_pythonize__(PyObject* klass, const std::string& name) {
          podio::detail::pythonizations::pythonize_subscript(klass, name);
  }
  /// Pythonizations for ROOT <= 6.40 before the switch to cppjit
  static void __cppyy_pythonize__(PyObject* klass, const std::string& name) {
        __cppjit_pythonize__(klass, name);
  }

  /// Print this collection to the passed stream
  void print(std::ostream& os=std::cout, bool flush=true) const final;

  /// Append a new object to the collection, and return this object.
  /// Does not invalidate iterators.
  MutableSenseWireHit create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableSenseWireHit create(Args&&... args);

  /// number of elements in the collection
  std::size_t size() const final;

  /// maximal number of elements in the collection
  std::size_t max_size() const final;

  /// Is the collection empty
  bool empty() const final;

  /// fully qualified type name
  const std::string_view getTypeName() const final { return typeName; }
  /// fully qualified type name of elements - with namespace
  const std::string_view getValueTypeName() const final { return valueTypeName; }
  /// fully qualified type name of stored POD elements - with namespace
  const std::string_view getDataTypeName() const final { return dataTypeName; }
  /// schema version
  podio::SchemaVersionT getSchemaVersion() const final;

  bool isSubsetCollection() const final {
    return m_isSubsetColl;
  }

  void setSubsetCollection(bool setSubset=true) final;

  /// Returns the const object of given index
  SenseWireHit operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableSenseWireHit operator[](std::size_t index);
  /// Returns the const object of given index
  SenseWireHit at(std::size_t index) const;
  /// Returns the object of given index
  MutableSenseWireHit at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableSenseWireHit& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const SenseWireHit& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (SenseWireHitObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
      );
    }
  }

  uint32_t getID() const final {
    return m_collectionID;
  }

  /// check if the collection has a valid ID
  bool hasID() const final {
    return getID() != static_cast<uint32_t>(podio::ObjectID::untracked) &&
        getID() != static_cast<uint32_t>(podio::ObjectID::invalid);
  }

  size_t getDatamodelRegistryIndex() const final;

  // support for the iterator protocol
  iterator begin() {
    return iterator(0, &m_storage.entries);
  }
  const_iterator begin() const {
    return const_iterator(0, &m_storage.entries);
  }
  const_iterator cbegin() const {
    return begin();
  }
  iterator end() {
    return iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator end() const {
    return const_iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator cend() const {
    return end();
  }
  // reverse iterators
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return rbegin();
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return rend();
  }


  std::vector<std::uint64_t> cellID(const size_t nElem = 0) const;
  std::vector<std::int32_t> type(const size_t nElem = 0) const;
  std::vector<std::int32_t> quality(const size_t nElem = 0) const;
  std::vector<float> time(const size_t nElem = 0) const;
  std::vector<float> eDep(const size_t nElem = 0) const;
  std::vector<float> eDepError(const size_t nElem = 0) const;
  std::vector<float> wireStereoAngle(const size_t nElem = 0) const;
  std::vector<float> wireAzimuthalAngle(const size_t nElem = 0) const;
  std::vector<edm4hep::Vector3d> position(const size_t nElem = 0) const;
  std::vector<double> positionAlongWireError(const size_t nElem = 0) const;
  std::vector<float> distanceToWire(const size_t nElem = 0) const;
  std::vector<float> distanceToWireError(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class SenseWireHitCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable SenseWireHitCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const SenseWireHitCollection& v);

template<typename... Args>
MutableSenseWireHit SenseWireHitCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new SenseWireHitObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  // Need to initialize the relation vectors manually for the {ObjectID, SenseWireHitData} constructor
  obj->m_nElectrons = new std::vector<std::uint16_t>();
  m_storage.createRelations(obj);
  return MutableSenseWireHit(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<SenseWireHitCollection::value_type> R>
SenseWireHitCollection::SenseWireHitCollection(std::from_range_t, R&& range)
  : SenseWireHitCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<SenseWireHitCollection::value_type> R>
SenseWireHitCollection SenseWireHitCollection::from(R&& range) {
  SenseWireHitCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHitCollection& collection);
#endif

/// @}
} // namespace extension


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitCollection::typeName;
constexpr std::string_view extension::SenseWireHitCollection::valueTypeName;
constexpr std::string_view extension::SenseWireHitCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitCollection::typeName;
constexpr std::string_view extension::SenseWireHitCollection::valueTypeName;
constexpr std::string_view extension::SenseWireHitCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLinkCollection_H
#define EXTENSION_SenseWireHitSimTrackerHitLinkCollection_H

// datamodel specific includes
#include "extension/SenseWireHitSimTrackerHitLink.h"
#include "extension/MutableSenseWireHitSimTrackerHitLink.h"
#include "extension/SenseWireHitSimTrackerHitLinkObj.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollectionData.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/detail/Pythonizations.h"
#include "podio/utilities/TypeHelpers.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <mutex>
#include <memory>
#include <cstddef>

namespace podio {
  struct RelationNames;
}

namespace extension {

/// \addtogroup extension_user_facing
/// @{

class SenseWireHitSimTrackerHitLinkCollectionIterator {
public:
  using value_type = SenseWireHitSimTrackerHitLink;
  using difference_type = ptrdiff_t;
  using reference = SenseWireHitSimTrackerHitLink;
  using pointer = SenseWireHitSimTrackerHitLink*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  SenseWireHitSimTrackerHitLinkCollectionIterator(size_t index, const SenseWireHitSimTrackerHitLinkObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr}), m_collection(collection) {}
  SenseWireHitSimTrackerHitLinkCollectionIterator() = default;

  SenseWireHitSimTrackerHitLinkCollectionIterator(const SenseWireHitSimTrackerHitLinkCollectionIterator&) = default;
  SenseWireHitSimTrackerHitLinkCollectionIterator(SenseWireHitSimTrackerHitLinkCollectionIterator&&) = default;
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator=(const SenseWireHitSimTrackerHitLinkCollectionIterator&) = default;
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator=(SenseWireHitSimTrackerHitLinkCollectionIterator&&) = default;
  ~SenseWireHitSimTrackerHitLinkCollectionIterator() = default;

  constexpr auto operator<=>(const SenseWireHitSimTrackerHitLinkCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const SenseWireHitSimTrackerHitLinkCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator++();
  SenseWireHitSimTrackerHitLinkCollectionIterator operator++(int);
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator--();
  SenseWireHitSimTrackerHitLinkCollectionIterator operator--(int);
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator+=(difference_type n);
  SenseWireHitSimTrackerHitLinkCollectionIterator operator+(difference_type n) const;
  friend SenseWireHitSimTrackerHitLinkCollectionIterator operator+(difference_type n, const SenseWireHitSimTrackerHitLinkCollectionIterator& it);
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator-=(difference_type n);
  SenseWireHitSimTrackerHitLinkCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const SenseWireHitSimTrackerHitLinkCollectionIterator& other) const;

private:
  size_t m_index{0};
  SenseWireHitSimTrackerHitLink m_object { podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr} };
  const SenseWireHitSimTrackerHitLinkObjPointerContainer* m_collection{nullptr};
};


class SenseWireHitSimTrackerHitLinkMutableCollectionIterator {
public:
  using value_type = SenseWireHitSimTrackerHitLink;
  using difference_type = ptrdiff_t;
  using reference = MutableSenseWireHitSimTrackerHitLink;
  using pointer = MutableSenseWireHitSimTrackerHitLink*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  SenseWireHitSimTrackerHitLinkMutableCollectionIterator(size_t index, const SenseWireHitSimTrackerHitLinkObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr}), m_collection(collection) {}
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator() = default;

  SenseWireHitSimTrackerHitLinkMutableCollectionIterator(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator&) = default;
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator(SenseWireHitSimTrackerHitLinkMutableCollectionIterator&&) = default;
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator=(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator&) = default;
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator=(SenseWireHitSimTrackerHitLinkMutableCollectionIterator&&) = default;
  ~SenseWireHitSimTrackerHitLinkMutableCollectionIterator() = default;

  constexpr auto operator<=>(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator++();
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator++(int);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator--();
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator--(int);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator+=(difference_type n);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator+(difference_type n) const;
  friend SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator+(difference_type n, const SenseWireHitSimTrackerHitLinkMutableCollectionIterator& it);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator-=(difference_type n);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableSenseWireHitSimTrackerHitLink m_object { podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr} };
  const SenseWireHitSimTrackerHitLinkObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class SenseWireHitSimTrackerHitLinkCollection : public podio::CollectionBase {
public:
  using value_type = SenseWireHitSimTrackerHitLink;
  using mutable_type = MutableSenseWireHitSimTrackerHitLink;
  using const_iterator = SenseWireHitSimTrackerHitLinkCollectionIterator;
  using iterator = SenseWireHitSimTrackerHitLinkMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  SenseWireHitSimTrackerHitLinkCollection() = default;
  SenseWireHitSimTrackerHitLinkCollection(SenseWireHitSimTrackerHitLinkCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  SenseWireHitSimTrackerHitLinkCollection(const SenseWireHitSimTrackerHitLinkCollection& ) = delete;
  SenseWireHitSimTrackerHitLinkCollection& operator=(const SenseWireHitSimTrackerHitLinkCollection& ) = delete;
  SenseWireHitSimTrackerHitLinkCollection(SenseWireHitSimTrackerHitLinkCollection&&) = default;
  SenseWireHitSimTrackerHitLinkCollection& operator=(SenseWireHitSimTrackerHitLinkCollection&&) = default;

  ~SenseWireHitSimTrackerHitLinkCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  SenseWireHitSimTrackerHitLinkCollection(std::from_range_t, R&& range);
#endif

  /// Construct a SenseWireHitSimTrackerHitLink collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static SenseWireHitSimTrackerHitLinkCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::SenseWireHitSimTrackerHitLinkCollection";
  constexpr static std::string_view valueTypeName = "extension::SenseWireHitSimTrackerHitLink";
  constexpr static std::string_view dataTypeName = "extension::SenseWireHitSimTrackerHitLinkData";

  /// Invalidates all iterators.
  void clear() final;

  /// Cppyy protocol to setup the pythonizations for this class. Not to be called directly.
  static void __cppjit_pythonize__(PyObject* klass, const std::string& name) {
          podio::detail::pythonizations::pythonize_subscript(klass, name);
  }
  /// Pythonizations for ROOT <= 6.40 before the switch to cppjit
  static void __cppyy_pythonize__(PyObject* klass, const std::string& name) {
        __cppjit_pythonize__(klass, name);
  }

  /// Print this collection to the passed stream
  void print(std::ostream& os=std::cout, bool flush=true) const final;

  /// Append a new object to the collection, and return this object.
  /// Does not invalidate iterators.
  MutableSenseWireHitSimTrackerHitLink create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableSenseWireHitSimTrackerHitLink create(Args&&... args);

  /// number of elements in the collection
  std::size_t size() const final;

  /// maximal number of elements in the collection
  std::size_t max_size() const final;

  /// Is the collection empty
  bool empty() const final;

  /// fully qualified type name
  const std::string_view getTypeName() const final { return typeName; }
  /// fully qualified type name of elements - with namespace
  const std::string_view getValueTypeName() const final { return valueTypeName; }
  /// fully qualified type name of stored POD elements - with namespace
  const std::string_view getDataTypeName() const final { return dataTypeName; }
  /// schema version
  podio::SchemaVersionT getSchemaVersion() const final;

  bool isSubsetCollection() const final {
    return m_isSubsetColl;
  }

  void setSubsetCollection(bool setSubset=true) final;

  /// Returns the const object of given index
  SenseWireHitSimTrackerHitLink operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableSenseWireHitSimTrackerHitLink operator[](std::size_t index);
  /// Returns the const object of given index
  SenseWireHitSimTrackerHitLink at(std::size_t index) const;
  /// Returns the object of given index
  MutableSenseWireHitSimTrackerHitLink at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableSenseWireHitSimTrackerHitLink& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const SenseWireHitSimTrackerHitLink& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (SenseWireHitSimTrackerHitLinkObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
      );
    }
  }

  uint32_t getID() const final {
    return m_collectionID;
  }

  /// check if the collection has a valid ID
  bool hasID() const final {
    return getID() != static_cast<uint32_t>(podio::ObjectID::untracked) &&
        getID() != static_cast<uint32_t>(podio::ObjectID::invalid);
  }

  size_t getDatamodelRegistryIndex() const final;

  // support for the iterator protocol
  iterator begin() {
    return iterator(0, &m_storage.entries);
  }
  const_iterator begin() const {
    return const_iterator(0, &m_storage.entries);
  }
  const_iterator cbegin() const {
    return begin();
  }
  iterator end() {
    return iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator end() const {
    return const_iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator cend() const {
    return end();
  }
  // reverse iterators
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return rbegin();
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return rend();
  }


  std::vector<float> weight(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class SenseWireHitSimTrackerHitLinkCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable SenseWireHitSimTrackerHitLinkCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const SenseWireHitSimTrackerHitLinkCollection& v);

template<typename... Args>
MutableSenseWireHitSimTrackerHitLink SenseWireHitSimTrackerHitLinkCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new SenseWireHitSimTrackerHitLinkObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  return MutableSenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<SenseWireHitSimTrackerHitLinkCollection::value_type> R>
SenseWireHitSimTrackerHitLinkCollection::SenseWireHitSimTrackerHitLinkCollection(std::from_range_t, R&& range)
  : SenseWireHitSimTrackerHitLinkCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<SenseWireHitSimTrackerHitLinkCollection::value_type> R>
SenseWireHitSimTrackerHitLinkCollection SenseWireHitSimTrackerHitLinkCollection::from(R&& range) {
  SenseWireHitSimTrackerHitLinkCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHitSimTrackerHitLinkCollection& collection);
#endif

/// @}
} // namespace extension


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::typeName;
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::valueTypeName;
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::typeName;
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::valueTypeName;
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_TrackCollection_H
#define EXTENSION_TrackCollection_H

// datamodel specific includes
#include "extension/Track.h"
#include "extension/MutableTrack.h"
#include "extension/TrackObj.h"
#include "extension/TrackCollectionData.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/detail/Pythonizations.h"
#include "podio/utilities/TypeHelpers.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <mutex>
#include <memory>
#include <cstddef>

namespace podio {
  struct RelationNames;
}

namespace extension {

/// \addtogroup extension_user_facing
/// @{

class TrackCollectionIterator {
public:
  using value_type = Track;
  using difference_type = ptrdiff_t;
  using reference = Track;
  using pointer = Track*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  TrackCollectionIterator(size_t index, const TrackObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<TrackObj>{nullptr}), m_collection(collection) {}
  TrackCollectionIterator() = default;

  TrackCollectionIterator(const TrackCollectionIterator&) = default;
  TrackCollectionIterator(TrackCollectionIterator&&) = default;
  TrackCollectionIterator& operator=(const TrackCollectionIterator&) = default;
  TrackCollectionIterator& operator=(TrackCollectionIterator&&) = default;
  ~TrackCollectionIterator() = default;

  constexpr auto operator<=>(const TrackCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const TrackCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  TrackCollectionIterator& operator++();
  TrackCollectionIterator operator++(int);
  TrackCollectionIterator& operator--();
  TrackCollectionIterator operator--(int);
  TrackCollectionIterator& operator+=(difference_type n);
  TrackCollectionIterator operator+(difference_type n) const;
  friend TrackCollectionIterator operator+(difference_type n, const TrackCollectionIterator& it);
  TrackCollectionIterator& operator-=(difference_type n);
  TrackCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const TrackCollectionIterator& other) const;

private:
  size_t m_index{0};
  Track m_object { podio::utils::MaybeSharedPtr<TrackObj>{nullptr} };
  const TrackObjPointerContainer* m_collection{nullptr};
};


class TrackMutableCollectionIterator {
public:
  using value_type = Track;
  using difference_type = ptrdiff_t;
  using reference = MutableTrack;
  using pointer = MutableTrack*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  TrackMutableCollectionIterator(size_t index, const TrackObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<TrackObj>{nullptr}), m_collection(collection) {}
  TrackMutableCollectionIterator() = default;

  TrackMutableCollectionIterator(const TrackMutableCollectionIterator&) = default;
  TrackMutableCollectionIterator(TrackMutableCollectionIterator&&) = default;
  TrackMutableCollectionIterator& operator=(const TrackMutableCollectionIterator&) = default;
  TrackMutableCollectionIterator& operator=(TrackMutableCollectionIterator&&) = default;
  ~TrackMutableCollectionIterator() = default;

  constexpr auto operator<=>(const TrackMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const TrackMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  TrackMutableCollectionIterator& operator++();
  TrackMutableCollectionIterator operator++(int);
  TrackMutableCollectionIterator& operator--();
  TrackMutableCollectionIterator operator--(int);
  TrackMutableCollectionIterator& operator+=(difference_type n);
  TrackMutableCollectionIterator operator+(difference_type n) const;
  friend TrackMutableCollectionIterator operator+(difference_type n, const TrackMutableCollectionIterator& it);
  TrackMutableCollectionIterator& operator-=(difference_type n);
  TrackMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const TrackMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableTrack m_object { podio::utils::MaybeSharedPtr<TrackObj>{nullptr} };
  const TrackObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class TrackCollection : public podio::CollectionBase {
public:
  using value_type = Track;
  using mutable_type = MutableTrack;
  using const_iterator = TrackCollectionIterator;
  using iterator = TrackMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  TrackCollection() = default;
  TrackCollection(TrackCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  TrackCollection(const TrackCollection& ) = delete;
  TrackCollection& operator=(const TrackCollection& ) = delete;
  TrackCollection(TrackCollection&&) = default;
  TrackCollection& operator=(TrackCollection&&) = default;

  ~TrackCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  TrackCollection(std::from_range_t, R&& range);
#endif

  /// Construct a Track collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static TrackCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::TrackCollection";
  constexpr static std::string_view valueTypeName = "extension::Track";
  constexpr static std::string_view dataTypeName = "extension::TrackData";

  /// Invalidates all iterators.
  void clear() final;

  /// Cppyy protocol to setup the pythonizations for this class. Not to be called directly.
  static void __cppjit_pythonize__(PyObject* klass, const std::string& name) {
          podio::detail::pythonizations::pythonize_subscript(klass, name);
  }
  /// Pythonizations for ROOT <= 6.40 before the switch to cppjit
  static void __cppyy_pythonize__(PyObject* klass, const std::string& name) {
        __cppjit_pythonize__(klass, name);
  }

  /// Print this collection to the passed stream
  void print(std::ostream& os=std::cout, bool flush=true) const final;

  /// Append a new object to the collection, and return this object.
  /// Does not invalidate iterators.
  MutableTrack create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableTrack create(Args&&... args);

  /// number of elements in the collection
  std::size_t size() const final;

  /// maximal number of elements in the collection
  std::size_t max_size() const final;

  /// Is the collection empty
  bool empty() const final;

  /// fully qualified type name
  const std::string_view getTypeName() const final { return typeName; }
  /// fully qualified type name of elements - with namespace
  const std::string_view getValueTypeName() const final { return valueTypeName; }
  /// fully qualified type name of stored POD elements - with namespace
  const std::string_view getDataTypeName() const final { return dataTypeName; }
  /// schema version
  podio::SchemaVersionT getSchemaVersion() const final;

  bool isSubsetCollection() const final {
    return m_isSubsetColl;
  }

  void setSubsetCollection(bool setSubset=true) final;

  /// Returns the const object of given index
  Track operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableTrack operator[](std::size_t index);
  /// Returns the const object of given index
  Track at(std::size_t index) const;
  /// Returns the object of given index
  MutableTrack at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableTrack& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const Track& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (TrackObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
      );
    }
  }

  uint32_t getID() const final {
    return m_collectionID;
  }

  /// check if the collection has a valid ID
  bool hasID() const final {
    return getID() != static_cast<uint32_t>(podio::ObjectID::untracked) &&
        getID() != static_cast<uint32_t>(podio::ObjectID::invalid);
  }

  size_t getDatamodelRegistryIndex() const final;

  // support for the iterator protocol
  iterator begin() {
    return iterator(0, &m_storage.entries);
  }
  const_iterator begin() const {
    return const_iterator(0, &m_storage.entries);
  }
  const_iterator cbegin() const {
    return begin();
  }
  iterator end() {
    return iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator end() const {
    return const_iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator cend() const {
    return end();
  }
  // reverse iterators
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return rbegin();
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return rend();
  }


  std::vector<std::int32_t> type(const size_t nElem = 0) const;
  std::vector<float> chi2(const size_t nElem = 0) const;
  std::vector<std::int32_t> ndf(const size_t nElem = 0) const;
  std::vector<float> dEdx(const size_t nElem = 0) const;
  std::vector<float> dEdxError(const size_t nElem = 0) const;
  std::vector<float> radiusOfInnermostHit(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class TrackCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable TrackCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const TrackCollection& v);

template<typename... Args>
MutableTrack TrackCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new TrackObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  // Need to initialize the relation vectors manually for the {ObjectID, TrackData} constructor
  obj->m_trackerHits = new std::vector<extension::TrackerHit>();
  obj->m_tracks = new std::vector<extension::Track>();
  obj->m_subdetectorHitNumbers = new std::vector<std::int32_t>();
  obj->m_trackStates = new std::vector<edm4hep::TrackState>();
  obj->m_dxQuantities = new std::vector<edm4hep::Quantity>();
  m_storage.createRelations(obj);
  return MutableTrack(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<TrackCollection::value_type> R>
TrackCollection::TrackCollection(std::from_range_t, R&& range)
  : TrackCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<TrackCollection::value_type> R>
TrackCollection TrackCollection::from(R&& range) {
  TrackCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const TrackCollection& collection);
#endif

/// @}
} // namespace extension


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::TrackCollection::typeName;
constexpr std::string_view extension::TrackCollection::valueTypeName;
constexpr std::string_view extension::TrackCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::TrackCollection::typeName;
constexpr std::string_view extension::TrackCollection::valueTypeName;
constexpr std::string_view extension::TrackCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef extensionDATAMODEL_DEFINITION_H
#define extensionDATAMODEL_DEFINITION_H

#include "podio/DatamodelRegistry.h"
#include "podio/SchemaEvolution.h"

/**
 * \defgroup extension_user_facing extension datamodel
 * \brief The user facing classes of the extension datamodel generated by podio.
 *
 * These are the classes (e.g. the immutable and `Mutable*` datatypes, their
 * `*Collection`s, components, interfaces and links) that make up the public
 * API of extension.
 */

/// \cond PODIO_INTERNAL
/**
 * \defgroup extension_podio_internal extension podio internals
 * \brief Implementation details of the extension datamodel generated by podio.
 *
 * These classes (e.g. the `*Obj`, `*Data`, `*CollectionData` and `*SIOBlock`
 * classes) are generated by podio for internal use only. They are not part of
 * the public API of extension and should not be used directly.
 */
/// \endcond

namespace extension::meta {
/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{
/**
 * The complete definition of the datamodel at generation time in JSON format.
 */
static constexpr auto extension__JSONDefinition = R"DATAMODELDEF({"options": {"getSyntax": true, "exposePODMembers": false, "includeSubfolder": "extension/", "setSyntax": true}, "schema_version": 1, "components": {}, "datatypes": {"extension::DriftChamberDigiLocal": {"Description": "Drift chamber digitized hit (before tracking) in local coordinates", "Author": "B. Francois, CERN", "Members": ["std::uint64_t cellID // ID of the wire that created this hit", "float distanceToWire // smeared distance of closest approach between the wire and the hit [mm]", "float zPositionAlongWire // smeared z position in the local wire coordinate system [mm]", "float time // time of the hit [ns].", "float eDep // energy deposited on the hit [GeV].", "float eDepError // error measured on eDep [GeV].", "std::uint32_t clusterCount // number of clusters associated to this hit"], "VectorMembers": [], "OneToOneRelations": [], "OneToManyRelations": [], "ExtraCode": {}, "MutableExtraCode": {}}, "extension::DriftChamberDigi": {"Description": "Drift chamber digitized hit (before tracking) in global coordinates. Assumes that the hits are radially in the middle of the cells", "Author": "B. Francois, CERN", "Members": ["std::uint64_t cellID // ID of the wire that created this hit", "edm4hep::Vector3d leftPosition // position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]", "edm4hep::Vector3d rightPosition // position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]", "float time // time of the hit [ns].", "float eDep // energy deposited on the hit [GeV].", "float eDepError // error measured on eDep [GeV].", "std::uint32_t clusterCount // number of clusters associated to this hit", "std::int32_t type // type of the raw data hit", "std::int32_t quality // quality bit flag of the hit", "edm4hep::Vector3d position[mm] // hit position"], "VectorMembers": [], "OneToOneRelations": [], "OneToManyRelations": [], "ExtraCode": {}, "MutableExtraCode": {}}, "extension::MCRecoDriftChamberDigiAssociation": {"Description": "Association between a DriftChamberDigi and the corresponding simulated hit", "Author": "B. Francois, CERN", "Members": ["float weight // weight of this association"], "OneToOneRelations": ["extension::DriftChamberDigi digi // reference to the digitized hit", "edm4hep::SimTrackerHit sim // reference to the simulated hit"], "VectorMembers": [], "OneToManyRelations": [], "ExtraCode": {}, "MutableExtraCode": {}}, "extension::DriftChamberDigiV2": {"Description": "Digitized hit (before tracking) for Drift Chamber v2 (requires data extension).", "Author": "A. Tolosa-Delgado, B. Francois, CERN", "Members": ["std::uint64_t cellID // ID of the sensor that created this hit", "std::int32_t type // type of the raw data hit", "std::int32_t quality // quality bit flag of the hit", "float time // time of the hit [ns]", "float eDep // energy deposited on the hit [GeV]", "float eDepError // error measured on eDep [GeV]", "edm4hep::Vector3d position // point on the sensitive wire (SW) which is closest to the simhit [mm]", "edm4hep::Vector3d directionSW // direction of SW", "float distanceToWire // distance hit-wire [mm]", "std::uint32_t nCluster // number of clusters associated to this hit"], "VectorMembers": ["std::uint16_t nElectrons // number of electrons for each cluster"], "OneToOneRelations": [], "OneToManyRelations": [], "ExtraCode": {}, "MutableExtraCode": {}}, "extension::MCRecoDriftChamberDigiV2Association": {"Description": "Association between a DriftChamberDigi and the corresponding simulated hit", "Author": "B. Francois, CERN", "Members": ["float weight // weight of this association"], "OneToOneRelations": ["extension::DriftChamberDigiV2 digi // reference to the digitized hit", "edm4hep::SimTrackerHit sim // reference to the simulated hit"], "VectorMembers": [], "OneToManyRelations": [], "ExtraCode": {}, "MutableExtraCode": {}}, "extension::Track": {"Description": "Reconstructed track", "Author": "EDM4hep authors", "Members": ["std::int32_t type // flagword that defines the type of track.Bits 16-31 are used internally", "float chi2 // Chi^2 of the track fit", "std::int32_t ndf // number of degrees of freedom of the track fit", "float dEdx // dEdx of the track", "float dEdxError // error of dEdx", "float radiusOfInnermostHit // radius of the innermost hit that has been used in the track fit"], "VectorMembers": ["std::int32_t subdetectorHitNumbers // number of hits in particular subdetectors", "edm4hep::TrackState trackStates // track states", "edm4hep::Quantity dxQuantities // different measurements of dx quantities"], "OneToManyRelations": ["extension::TrackerHit trackerHits // hits that have been used to create this track", "extension::Track tracks // tracks (segments) that have been combined to create this track"], "OneToOneRelations": [], "ExtraCode": {}, "MutableExtraCode": {}}, "extension::SenseWireHit": {"Description": "Sense wire hit, before the hit position degeneracy is resolved. The circle representing possible positions is parametrized with its center, radius and normal vector (given by the wire direction).", "Author": "EDM4hep authors", "Members": ["std::uint64_t cellID // ID of the sensor that created this hit", "std::int32_t type // type of the raw data hit", "std::int32_t quality // quality bit flag of the hit", "float time[ns] // time of the hit", "float eDep[GeV] // energy deposited by the hit", "float eDepError[GeV] // error on eDep", "float wireStereoAngle // angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction", "float wireAzimuthalAngle // azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction", "edm4hep::Vector3d position[mm] // point on the sense wire which is closest to the hit (center of the circle)", "double positionAlongWireError[mm] // error on the hit position along the wire direction", "float distanceToWire[mm] // distance between the hit and the wire (radius of the circle)", "float distanceToWireError[mm] // error on distanceToWire"], "VectorMembers": ["std::uint16_t nElectrons // number of electrons for each cluster (number of clusters = vector size)"], "ExtraCode": {"declaration": " /// Return the number of clusters associated to the hit\n auto getNClusters() const { return getNElectrons().size(); }\n "}, "OneToOneRelations": [], "OneToManyRelations": [], "MutableExtraCode": {}}, "extension::SenseWireHitSimTrackerHitLink": {"Description": "Link between a SenseWireHit and a SimTrackerHit", "Author": "EDM4hep authors", "Members": ["float weight // weight of this link"], "OneToOneRelations": ["extension::SenseWireHit from // reference to the SenseWireHit", "edm4hep::SimTrackerHit to // reference to the SimTrackerHit"], "VectorMembers": [], "OneToManyRelations": [], "ExtraCode": {}, "MutableExtraCode": {}}}, "interfaces": {"extension::TrackerHit": {"Description": "Tracker hit interface class", "Author": "Thomas Madlener, DESY", "Members": ["std::uint64_t cellID // ID of the sensor that created this hit", "std::int32_t type // type of the raw data hit", "std::int32_t quality // quality bit flag of the hit", "float time[ns] // time of the hit", "float eDep[GeV] // energy deposited on the hit", "float eDepError[GeV] // error measured on eDep", "edm4hep::Vector3d position[mm] // hit position"], "Types": ["edm4hep::TrackerHitPlane", "edm4hep::TrackerHit3D", "extension::SenseWireHit", "extension::DriftChamberDigi"]}, "extension::TrackerHit_dev": {"Description": "Tracker hit interface class to test the mixing of extension and edm4hep data types", "Author": "Thomas Madlener, DESY", "Members": ["std::uint64_t cellID // ID of the sensor that created this hit", "std::int32_t type // type of the raw data hit", "std::int32_t quality // quality bit flag of the hit", "float time[ns] // time of the hit", "float eDep[GeV] // energy deposited on the hit", "float eDepError[GeV] // error measured on eDep", "edm4hep::Vector3d position[mm] // hit position"], "Types": ["edm4hep::TrackerHitPlane", "edm4hep::TrackerHit3D", "extension::SenseWireHit"]}}, "links": {}})DATAMODELDEF";


/**
 * The names of all relations and vector members for all datatypes
 */
inline podio::RelationNameMapping extension__getRelationNames() {
  using namespace std::string_view_literals;
  return {
    {"extension::DriftChamberDigiLocal"sv,
      {  },
      {  },
    },
    {"extension::DriftChamberDigi"sv,
      {  },
      {  },
    },
    {"extension::MCRecoDriftChamberDigiAssociation"sv,
      { "digi"sv, "sim"sv },
      {  },
    },
    {"extension::DriftChamberDigiV2"sv,
      {  },
      { "nElectrons"sv },
    },
    {"extension::MCRecoDriftChamberDigiV2Association"sv,
      { "digi"sv, "sim"sv },
      {  },
    },
    {"extension::Track"sv,
      { "trackerHits"sv, "tracks"sv },
      { "subdetectorHitNumbers"sv, "trackStates"sv, "dxQuantities"sv },
    },
    {"extension::SenseWireHit"sv,
      {  },
      { "nElectrons"sv },
    },
    {"extension::SenseWireHitSimTrackerHitLink"sv,
      { "from"sv, "to"sv },
      {  },
    },
  };
}
/// @}
/// \endcond

/**
 * The schema version at generation time
 */
static constexpr podio::SchemaVersionT schemaVersion = 1;

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{
/**
 * The helper class that takes care of registering the datamodel definition to
 * the DatamodelRegistry and to provide the index in that registry.
 *
 * Implemented as a singleton mainly to ensure only a single registration of
 * each datamodel, during the constructor
 */
class DatamodelRegistryIndex {
public:
  static size_t value() {
    static const auto relationNames = extension__getRelationNames();
    static auto index =
      DatamodelRegistryIndex(podio::DatamodelRegistry::mutInstance().registerDatamodel(
                               "extension",
                               extension__JSONDefinition,
                               relationNames
                             ));
    return index.m_value;
  }
private:
  DatamodelRegistryIndex(size_t v) : m_value(v) {}
  size_t m_value{podio::DatamodelRegistry::NoDefinitionAvailable};
};


namespace static_registration {
  // The usual trick via an IIFE and a const variable that we assign to, to
  // ensure that we populate this before everything starts
  inline  bool ensureRegistration() {
    const static auto reg = []() {
      return extension::meta::DatamodelRegistryIndex::value() != podio::DatamodelRegistry::NoDefinitionAvailable;
    }();
    return reg;
  }

  const auto registrationEnsured = ensureRegistration();
}
/// @}
/// \endcond

} // namespace extension::meta

#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigi_H
#define EXTENSION_DriftChamberDigi_H

#include "extension/DriftChamberDigiObj.h"

#include "edm4hep/Vector3d.h"
#include <cstdint>

#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class DriftChamberDigiCollection;
class TrackerHit;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::DriftChamberDigi& obj);
};

namespace extension {

class MutableDriftChamberDigi;
class DriftChamberDigiCollection;
class DriftChamberDigiCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class DriftChamberDigi
 *  Drift chamber digitized hit (before tracking) in global coordinates. Assumes that the hits are radially in the middle of the cells
 *  @author: B. Francois, CERN
 */
class DriftChamberDigi {

  friend class MutableDriftChamberDigi;
  friend class DriftChamberDigiCollection;
  friend class extension::DriftChamberDigiCollectionData;
  friend class DriftChamberDigiCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const DriftChamberDigi & obj);

public:
  using mutable_type = MutableDriftChamberDigi;
  using collection_type = DriftChamberDigiCollection;

  /// default constructor
  DriftChamberDigi() = default;

  /// Constructor initializing all members
  DriftChamberDigi(const std::uint64_t cellID, const edm4hep::Vector3d& leftPosition, const edm4hep::Vector3d& rightPosition, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount, const std::int32_t type, const std::int32_t quality, const edm4hep::Vector3d& position);

  /// copy constructor
  DriftChamberDigi(const DriftChamberDigi& other) noexcept = default;

  /// move constructor
  DriftChamberDigi(DriftChamberDigi&& other) noexcept = default;

  /// copy-assignment operator
  DriftChamberDigi& operator=(DriftChamberDigi other) &; // Rebind this to other's internal object
  DriftChamberDigi& operator=(DriftChamberDigi other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigi clone(bool cloneRelations=true) const;

  /// destructor
  ~DriftChamberDigi() = default;

  /// converting constructor from mutable object
  DriftChamberDigi(const MutableDriftChamberDigi& other);

  static DriftChamberDigi makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::DriftChamberDigi";

  /// Access the ID of the wire that created this hit
  std::uint64_t getCellID() const;

  /// Access the position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]
  const edm4hep::Vector3d& getLeftPosition() const;

  /// Access the position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]
  const edm4hep::Vector3d& getRightPosition() const;

  /// Access the time of the hit [ns].
  float getTime() const;

  /// Access the energy deposited on the hit [GeV].
  float getEDep() const;

  /// Access the error measured on eDep [GeV].
  float getEDepError() const;

  /// Access the number of clusters associated to this hit
  std::uint32_t getClusterCount() const;

  /// Access the type of the raw data hit
  std::int32_t getType() const;

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const;

  /// Access the hit position [mm]
  const edm4hep::Vector3d& getPosition() const;





  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiObj>{nullptr}; }

  constexpr bool operator==(const DriftChamberDigi& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableDriftChamberDigi& other) const;

  constexpr bool operator!=(const DriftChamberDigi& other) const { return !(*this == other); }
  bool operator!=(const MutableDriftChamberDigi& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const DriftChamberDigi& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<DriftChamberDigi>;

  friend void swap(DriftChamberDigi& a, DriftChamberDigi& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiObj
  explicit DriftChamberDigi(podio::utils::MaybeSharedPtr<DriftChamberDigiObj> obj);
  DriftChamberDigi(DriftChamberDigiObj* obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiObj> m_obj{new DriftChamberDigiObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const DriftChamberDigi& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigi& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::DriftChamberDigi> {
  std::size_t operator()(const extension::DriftChamberDigi& obj) const {
    return std::hash<extension::DriftChamberDigiObj*>{}(obj.m_obj.get());
  }
};


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigi::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigi::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigi_CollectionData_H
#define EXTENSION_DriftChamberDigi_CollectionData_H

// datamodel specific includes
#include "extension/DriftChamberDigiData.h"
#include "extension/DriftChamberDigiObj.h"
#include "extension/DriftChamberDigi.h"

// schema evolution specific includes

// podio specific includes
#include "podio/CollectionBuffers.h"
#include "podio/ICollectionProvider.h"

#include <deque>
#include <memory>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


using DriftChamberDigiObjPointerContainer = std::deque<DriftChamberDigiObj*>;
using DriftChamberDigiDataContainer = std::vector<DriftChamberDigiData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class DriftChamberDigiCollectionData {
public:
  /**
   * The Objs of this collection
   */
  DriftChamberDigiObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  DriftChamberDigiCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  DriftChamberDigiCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  DriftChamberDigiCollectionData(const DriftChamberDigiCollectionData&) = delete;
  DriftChamberDigiCollectionData& operator=(const DriftChamberDigiCollectionData&) = delete;
  DriftChamberDigiCollectionData(DriftChamberDigiCollectionData&& other) = default;
  DriftChamberDigiCollectionData& operator=(DriftChamberDigiCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~DriftChamberDigiCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();


  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations

  // members to handle vector members

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<DriftChamberDigiDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiDATA_H
#define EXTENSION_DriftChamberDigiDATA_H

#include "edm4hep/Vector3d.h"
#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class DriftChamberDigiData
 *  Drift chamber digitized hit (before tracking) in global coordinates. Assumes that the hits are radially in the middle of the cells
 *  @author: B. Francois, CERN
 */
class DriftChamberDigiData {
public:
  std::uint64_t cellID{}; ///< ID of the wire that created this hit
  ::edm4hep::Vector3d leftPosition{}; ///< position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]
  ::edm4hep::Vector3d rightPosition{}; ///< position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]
  float time{}; ///< time of the hit [ns].
  float eDep{}; ///< energy deposited on the hit [GeV].
  float eDepError{}; ///< error measured on eDep [GeV].
  std::uint32_t clusterCount{}; ///< number of clusters associated to this hit
  std::int32_t type{}; ///< type of the raw data hit
  std::int32_t quality{}; ///< quality bit flag of the hit
  ::edm4hep::Vector3d position{}; ///< hit position [mm]

};


namespace v1 {
using DriftChamberDigiData = extension::DriftChamberDigiData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiLocal_H
#define EXTENSION_DriftChamberDigiLocal_H

#include "extension/DriftChamberDigiLocalObj.h"

#include <cstdint>

#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class DriftChamberDigiLocalCollection;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::DriftChamberDigiLocal& obj);
};

namespace extension {

class MutableDriftChamberDigiLocal;
class DriftChamberDigiLocalCollection;
class DriftChamberDigiLocalCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class DriftChamberDigiLocal
 *  Drift chamber digitized hit (before tracking) in local coordinates
 *  @author: B. Francois, CERN
 */
class DriftChamberDigiLocal {

  friend class MutableDriftChamberDigiLocal;
  friend class DriftChamberDigiLocalCollection;
  friend class extension::DriftChamberDigiLocalCollectionData;
  friend class DriftChamberDigiLocalCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const DriftChamberDigiLocal & obj);

public:
  using mutable_type = MutableDriftChamberDigiLocal;
  using collection_type = DriftChamberDigiLocalCollection;

  /// default constructor
  DriftChamberDigiLocal() = default;

  /// Constructor initializing all members
  DriftChamberDigiLocal(const std::uint64_t cellID, const float distanceToWire, const float zPositionAlongWire, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount);

  /// copy constructor
  DriftChamberDigiLocal(const DriftChamberDigiLocal& other) noexcept = default;

  /// move constructor
  DriftChamberDigiLocal(DriftChamberDigiLocal&& other) noexcept = default;

  /// copy-assignment operator
  DriftChamberDigiLocal& operator=(DriftChamberDigiLocal other) &; // Rebind this to other's internal object
  DriftChamberDigiLocal& operator=(DriftChamberDigiLocal other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigiLocal clone(bool cloneRelations=true) const;

  /// destructor
  ~DriftChamberDigiLocal() = default;

  /// converting constructor from mutable object
  DriftChamberDigiLocal(const MutableDriftChamberDigiLocal& other);

  static DriftChamberDigiLocal makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::DriftChamberDigiLocal";

  /// Access the ID of the wire that created this hit
  std::uint64_t getCellID() const;

  /// Access the smeared distance of closest approach between the wire and the hit [mm]
  float getDistanceToWire() const;

  /// Access the smeared z position in the local wire coordinate system [mm]
  float getZPositionAlongWire() const;

  /// Access the time of the hit [ns].
  float getTime() const;

  /// Access the energy deposited on the hit [GeV].
  float getEDep() const;

  /// Access the error measured on eDep [GeV].
  float getEDepError() const;

  /// Access the number of clusters associated to this hit
  std::uint32_t getClusterCount() const;





  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiLocalObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>{nullptr}; }

  constexpr bool operator==(const DriftChamberDigiLocal& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableDriftChamberDigiLocal& other) const;

  constexpr bool operator!=(const DriftChamberDigiLocal& other) const { return !(*this == other); }
  bool operator!=(const MutableDriftChamberDigiLocal& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const DriftChamberDigiLocal& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<DriftChamberDigiLocal>;

  friend void swap(DriftChamberDigiLocal& a, DriftChamberDigiLocal& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiLocalObj
  explicit DriftChamberDigiLocal(podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> obj);
  DriftChamberDigiLocal(DriftChamberDigiLocalObj* obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> m_obj{new DriftChamberDigiLocalObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const DriftChamberDigiLocal& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiLocal& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::DriftChamberDigiLocal> {
  std::size_t operator()(const extension::DriftChamberDigiLocal& obj) const {
    return std::hash<extension::DriftChamberDigiLocalObj*>{}(obj.m_obj.get());
  }
};


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiLocal::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiLocal::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiLocal_CollectionData_H
#define EXTENSION_DriftChamberDigiLocal_CollectionData_H

// datamodel specific includes
#include "extension/DriftChamberDigiLocalData.h"
#include "extension/DriftChamberDigiLocalObj.h"
#include "extension/DriftChamberDigiLocal.h"

// schema evolution specific includes

// podio specific includes
#include "podio/CollectionBuffers.h"
#include "podio/ICollectionProvider.h"

#include <deque>
#include <memory>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


using DriftChamberDigiLocalObjPointerContainer = std::deque<DriftChamberDigiLocalObj*>;
using DriftChamberDigiLocalDataContainer = std::vector<DriftChamberDigiLocalData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class DriftChamberDigiLocalCollectionData {
public:
  /**
   * The Objs of this collection
   */
  DriftChamberDigiLocalObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  DriftChamberDigiLocalCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  DriftChamberDigiLocalCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  DriftChamberDigiLocalCollectionData(const DriftChamberDigiLocalCollectionData&) = delete;
  DriftChamberDigiLocalCollectionData& operator=(const DriftChamberDigiLocalCollectionData&) = delete;
  DriftChamberDigiLocalCollectionData(DriftChamberDigiLocalCollectionData&& other) = default;
  DriftChamberDigiLocalCollectionData& operator=(DriftChamberDigiLocalCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~DriftChamberDigiLocalCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();


  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations

  // members to handle vector members

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<DriftChamberDigiLocalDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiLocalDATA_H
#define EXTENSION_DriftChamberDigiLocalDATA_H

#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class DriftChamberDigiLocalData
 *  Drift chamber digitized hit (before tracking) in local coordinates
 *  @author: B. Francois, CERN
 */
class DriftChamberDigiLocalData {
public:
  std::uint64_t cellID{}; ///< ID of the wire that created this hit
  float distanceToWire{}; ///< smeared distance of closest approach between the wire and the hit [mm]
  float zPositionAlongWire{}; ///< smeared z position in the local wire coordinate system [mm]
  float time{}; ///< time of the hit [ns].
  float eDep{}; ///< energy deposited on the hit [GeV].
  float eDepError{}; ///< error measured on eDep [GeV].
  std::uint32_t clusterCount{}; ///< number of clusters associated to this hit

};


namespace v1 {
using DriftChamberDigiLocalData = extension::DriftChamberDigiLocalData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiLocalOBJ_H
#define EXTENSION_DriftChamberDigiLocalOBJ_H

// data model specific includes
#include "extension/DriftChamberDigiLocalData.h"

#include "podio/ObjectID.h"



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class DriftChamberDigiLocal;

class DriftChamberDigiLocalObj {
public:
  /// constructor
  DriftChamberDigiLocalObj();
  /// copy constructor (does a deep-copy of relation containers)
  DriftChamberDigiLocalObj(const DriftChamberDigiLocalObj&);
  /// constructor from ObjectID and DriftChamberDigiLocalData
  /// does not initialize the internal relation containers
  DriftChamberDigiLocalObj(const podio::ObjectID& id, const DriftChamberDigiLocalData& data);
  /// No assignment operator
  DriftChamberDigiLocalObj& operator=(const DriftChamberDigiLocalObj&) = delete;
  virtual ~DriftChamberDigiLocalObj() = default;

public:
  podio::ObjectID id{};
  DriftChamberDigiLocalData data;
};

/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiOBJ_H
#define EXTENSION_DriftChamberDigiOBJ_H

// data model specific includes
#include "extension/DriftChamberDigiData.h"

#include "podio/ObjectID.h"



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class DriftChamberDigi;

class DriftChamberDigiObj {
public:
  /// constructor
  DriftChamberDigiObj();
  /// copy constructor (does a deep-copy of relation containers)
  DriftChamberDigiObj(const DriftChamberDigiObj&);
  /// constructor from ObjectID and DriftChamberDigiData
  /// does not initialize the internal relation containers
  DriftChamberDigiObj(const podio::ObjectID& id, const DriftChamberDigiData& data);
  /// No assignment operator
  DriftChamberDigiObj& operator=(const DriftChamberDigiObj&) = delete;
  virtual ~DriftChamberDigiObj() = default;

public:
  podio::ObjectID id{};
  DriftChamberDigiData data;
};

/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2_H
#define EXTENSION_DriftChamberDigiV2_H

#include "extension/DriftChamberDigiV2Obj.h"

#include "edm4hep/Vector3d.h"
#include "podio/RelationRange.h"
#include <cstdint>
#include <vector>

#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class DriftChamberDigiV2Collection;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::DriftChamberDigiV2& obj);
};

namespace extension {

class MutableDriftChamberDigiV2;
class DriftChamberDigiV2Collection;
class DriftChamberDigiV2CollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class DriftChamberDigiV2
 *  Digitized hit (before tracking) for Drift Chamber v2 (requires data extension).
 *  @author: A. Tolosa-Delgado, B. Francois, CERN
 */
class DriftChamberDigiV2 {

  friend class MutableDriftChamberDigiV2;
  friend class DriftChamberDigiV2Collection;
  friend class extension::DriftChamberDigiV2CollectionData;
  friend class DriftChamberDigiV2CollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const DriftChamberDigiV2 & obj);

public:
  using mutable_type = MutableDriftChamberDigiV2;
  using collection_type = DriftChamberDigiV2Collection;

  /// default constructor
  DriftChamberDigiV2() = default;

  /// Constructor initializing all members
  DriftChamberDigiV2(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const edm4hep::Vector3d& position, const edm4hep::Vector3d& directionSW, const float distanceToWire, const std::uint32_t nCluster);

  /// copy constructor
  DriftChamberDigiV2(const DriftChamberDigiV2& other) noexcept = default;

  /// move constructor
  DriftChamberDigiV2(DriftChamberDigiV2&& other) noexcept = default;

  /// copy-assignment operator
  DriftChamberDigiV2& operator=(DriftChamberDigiV2 other) &; // Rebind this to other's internal object
  DriftChamberDigiV2& operator=(DriftChamberDigiV2 other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigiV2 clone(bool cloneRelations=true) const;

  /// destructor
  ~DriftChamberDigiV2() = default;

  /// converting constructor from mutable object
  DriftChamberDigiV2(const MutableDriftChamberDigiV2& other);

  static DriftChamberDigiV2 makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::DriftChamberDigiV2";

  /// Access the ID of the sensor that created this hit
  std::uint64_t getCellID() const;

  /// Access the type of the raw data hit
  std::int32_t getType() const;

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const;

  /// Access the time of the hit [ns]
  float getTime() const;

  /// Access the energy deposited on the hit [GeV]
  float getEDep() const;

  /// Access the error measured on eDep [GeV]
  float getEDepError() const;

  /// Access the point on the sensitive wire (SW) which is closest to the simhit [mm]
  const edm4hep::Vector3d& getPosition() const;

  /// Access the direction of SW
  const edm4hep::Vector3d& getDirectionSW() const;

  /// Access the distance hit-wire [mm]
  float getDistanceToWire() const;

  /// Access the number of clusters associated to this hit
  std::uint32_t getNCluster() const;



  std::size_t nElectrons_size() const;
  std::uint16_t getNElectrons(std::size_t) const;
  std::vector<std::uint16_t>::const_iterator nElectrons_begin() const;
  std::vector<std::uint16_t>::const_iterator nElectrons_end() const;
  podio::RelationRange<std::uint16_t> getNElectrons() const;


  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiV2Obj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr}; }

  constexpr bool operator==(const DriftChamberDigiV2& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableDriftChamberDigiV2& other) const;

  constexpr bool operator!=(const DriftChamberDigiV2& other) const { return !(*this == other); }
  bool operator!=(const MutableDriftChamberDigiV2& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const DriftChamberDigiV2& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<DriftChamberDigiV2>;

  friend void swap(DriftChamberDigiV2& a, DriftChamberDigiV2& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiV2Obj
  explicit DriftChamberDigiV2(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj> obj);
  DriftChamberDigiV2(DriftChamberDigiV2Obj* obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj> m_obj{new DriftChamberDigiV2Obj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const DriftChamberDigiV2& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiV2& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::DriftChamberDigiV2> {
  std::size_t operator()(const extension::DriftChamberDigiV2& obj) const {
    return std::hash<extension::DriftChamberDigiV2Obj*>{}(obj.m_obj.get());
  }
};


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiV2::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiV2::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2_CollectionData_H
#define EXTENSION_DriftChamberDigiV2_CollectionData_H

// datamodel specific includes
#include "extension/DriftChamberDigiV2Data.h"
#include "extension/DriftChamberDigiV2Obj.h"
#include "extension/DriftChamberDigiV2.h"

// schema evolution specific includes

// podio specific includes
#include "podio/CollectionBuffers.h"
#include "podio/ICollectionProvider.h"

#include <deque>
#include <memory>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


using DriftChamberDigiV2ObjPointerContainer = std::deque<DriftChamberDigiV2Obj*>;
using DriftChamberDigiV2DataContainer = std::vector<DriftChamberDigiV2Data>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class DriftChamberDigiV2CollectionData {
public:
  /**
   * The Objs of this collection
   */
  DriftChamberDigiV2ObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  DriftChamberDigiV2CollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  DriftChamberDigiV2CollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  DriftChamberDigiV2CollectionData(const DriftChamberDigiV2CollectionData&) = delete;
  DriftChamberDigiV2CollectionData& operator=(const DriftChamberDigiV2CollectionData&) = delete;
  DriftChamberDigiV2CollectionData(DriftChamberDigiV2CollectionData&& other) = default;
  DriftChamberDigiV2CollectionData& operator=(DriftChamberDigiV2CollectionData&& other) = default;

  /**
   * Destructor
   */
  ~DriftChamberDigiV2CollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();

  void createRelations(DriftChamberDigiV2Obj* obj);

  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations

  // members to handle vector members
  podio::UVecPtr<std::uint16_t> m_vec_nElectrons{nullptr}; /// combined vector of all objects in collection
  std::vector<podio::UVecPtr<std::uint16_t>> m_vecs_nElectrons{}; /// pointers to individual member vectors

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<DriftChamberDigiV2DataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2DATA_H
#define EXTENSION_DriftChamberDigiV2DATA_H

#include "edm4hep/Vector3d.h"
#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class DriftChamberDigiV2Data
 *  Digitized hit (before tracking) for Drift Chamber v2 (requires data extension).
 *  @author: A. Tolosa-Delgado, B. Francois, CERN
 */
class DriftChamberDigiV2Data {
public:
  std::uint64_t cellID{}; ///< ID of the sensor that created this hit
  std::int32_t type{}; ///< type of the raw data hit
  std::int32_t quality{}; ///< quality bit flag of the hit
  float time{}; ///< time of the hit [ns]
  float eDep{}; ///< energy deposited on the hit [GeV]
  float eDepError{}; ///< error measured on eDep [GeV]
  ::edm4hep::Vector3d position{}; ///< point on the sensitive wire (SW) which is closest to the simhit [mm]
  ::edm4hep::Vector3d directionSW{}; ///< direction of SW
  float distanceToWire{}; ///< distance hit-wire [mm]
  std::uint32_t nCluster{}; ///< number of clusters associated to this hit

  unsigned int nElectrons_begin{};
  unsigned int nElectrons_end{};
};


namespace v1 {
using DriftChamberDigiV2Data = extension::DriftChamberDigiV2Data;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2OBJ_H
#define EXTENSION_DriftChamberDigiV2OBJ_H

// data model specific includes
#include "extension/DriftChamberDigiV2Data.h"
#include <vector>

#include "podio/ObjectID.h"
#include <vector>



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class DriftChamberDigiV2;

class DriftChamberDigiV2Obj {
public:
  /// constructor
  DriftChamberDigiV2Obj();
  /// copy constructor (does a deep-copy of relation containers)
  DriftChamberDigiV2Obj(const DriftChamberDigiV2Obj&);
  /// constructor from ObjectID and DriftChamberDigiV2Data
  /// does not initialize the internal relation containers
  DriftChamberDigiV2Obj(const podio::ObjectID& id, const DriftChamberDigiV2Data& data);
  /// No assignment operator
  DriftChamberDigiV2Obj& operator=(const DriftChamberDigiV2Obj&) = delete;
  virtual ~DriftChamberDigiV2Obj();

public:
  podio::ObjectID id{};
  DriftChamberDigiV2Data data;
  std::vector<std::uint16_t>* m_nElectrons{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociation_H
#define EXTENSION_MCRecoDriftChamberDigiAssociation_H

#include "extension/MCRecoDriftChamberDigiAssociationObj.h"


#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class MCRecoDriftChamberDigiAssociationCollection;
class DriftChamberDigi;
class MutableDriftChamberDigi;
}
namespace edm4hep {
class SimTrackerHit;
class MutableSimTrackerHit;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::MCRecoDriftChamberDigiAssociation& obj);
};

namespace extension {

class MutableMCRecoDriftChamberDigiAssociation;
class MCRecoDriftChamberDigiAssociationCollection;
class MCRecoDriftChamberDigiAssociationCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class MCRecoDriftChamberDigiAssociation
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MCRecoDriftChamberDigiAssociation {

  friend class MutableMCRecoDriftChamberDigiAssociation;
  friend class MCRecoDriftChamberDigiAssociationCollection;
  friend class extension::MCRecoDriftChamberDigiAssociationCollectionData;
  friend class MCRecoDriftChamberDigiAssociationCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const MCRecoDriftChamberDigiAssociation & obj);

public:
  using mutable_type = MutableMCRecoDriftChamberDigiAssociation;
  using collection_type = MCRecoDriftChamberDigiAssociationCollection;

  /// default constructor
  MCRecoDriftChamberDigiAssociation() = default;

  /// Constructor initializing all members
  MCRecoDriftChamberDigiAssociation(const float weight);

  /// copy constructor
  MCRecoDriftChamberDigiAssociation(const MCRecoDriftChamberDigiAssociation& other) noexcept = default;

  /// move constructor
  MCRecoDriftChamberDigiAssociation(MCRecoDriftChamberDigiAssociation&& other) noexcept = default;

  /// copy-assignment operator
  MCRecoDriftChamberDigiAssociation& operator=(MCRecoDriftChamberDigiAssociation other) &; // Rebind this to other's internal object
  MCRecoDriftChamberDigiAssociation& operator=(MCRecoDriftChamberDigiAssociation other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableMCRecoDriftChamberDigiAssociation clone(bool cloneRelations=true) const;

  /// destructor
  ~MCRecoDriftChamberDigiAssociation() = default;

  /// converting constructor from mutable object
  MCRecoDriftChamberDigiAssociation(const MutableMCRecoDriftChamberDigiAssociation& other);

  static MCRecoDriftChamberDigiAssociation makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::MCRecoDriftChamberDigiAssociation";

  /// Access the weight of this association
  float getWeight() const;


  /// Access the reference to the digitized hit
  const extension::DriftChamberDigi getDigi() const;
  /// Access the reference to the simulated hit
  const edm4hep::SimTrackerHit getSim() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from MCRecoDriftChamberDigiAssociationObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr}; }

  constexpr bool operator==(const MCRecoDriftChamberDigiAssociation& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableMCRecoDriftChamberDigiAssociation& other) const;

  constexpr bool operator!=(const MCRecoDriftChamberDigiAssociation& other) const { return !(*this == other); }
  bool operator!=(const MutableMCRecoDriftChamberDigiAssociation& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MCRecoDriftChamberDigiAssociation& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MCRecoDriftChamberDigiAssociation>;

  friend void swap(MCRecoDriftChamberDigiAssociation& a, MCRecoDriftChamberDigiAssociation& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing MCRecoDriftChamberDigiAssociationObj
  explicit MCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> obj);
  MCRecoDriftChamberDigiAssociation(MCRecoDriftChamberDigiAssociationObj* obj);

  podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> m_obj{new MCRecoDriftChamberDigiAssociationObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiAssociation& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiAssociation& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MCRecoDriftChamberDigiAssociation> {
  std::size_t operator()(const extension::MCRecoDriftChamberDigiAssociation& obj) const {
    return std::hash<extension::MCRecoDriftChamberDigiAssociationObj*>{}(obj.m_obj.get());
  }
};


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociation::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociation::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociation_CollectionData_H
#define EXTENSION_MCRecoDriftChamberDigiAssociation_CollectionData_H

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiAssociationData.h"
#include "extension/MCRecoDriftChamberDigiAssociationObj.h"
#include "extension/MCRecoDriftChamberDigiAssociation.h"
#include "extension/DriftChamberDigi.h"
#include "edm4hep/SimTrackerHit.h"

// schema evolution specific includes

// podio specific includes
#include "podio/CollectionBuffers.h"
#include "podio/ICollectionProvider.h"

#include <deque>
#include <memory>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


using MCRecoDriftChamberDigiAssociationObjPointerContainer = std::deque<MCRecoDriftChamberDigiAssociationObj*>;
using MCRecoDriftChamberDigiAssociationDataContainer = std::vector<MCRecoDriftChamberDigiAssociationData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class MCRecoDriftChamberDigiAssociationCollectionData {
public:
  /**
   * The Objs of this collection
   */
  MCRecoDriftChamberDigiAssociationObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  MCRecoDriftChamberDigiAssociationCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  MCRecoDriftChamberDigiAssociationCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  MCRecoDriftChamberDigiAssociationCollectionData(const MCRecoDriftChamberDigiAssociationCollectionData&) = delete;
  MCRecoDriftChamberDigiAssociationCollectionData& operator=(const MCRecoDriftChamberDigiAssociationCollectionData&) = delete;
  MCRecoDriftChamberDigiAssociationCollectionData(MCRecoDriftChamberDigiAssociationCollectionData&& other) = default;
  MCRecoDriftChamberDigiAssociationCollectionData& operator=(MCRecoDriftChamberDigiAssociationCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~MCRecoDriftChamberDigiAssociationCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();


  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations
  podio::UVecPtr<extension::DriftChamberDigi> m_rel_digi{std::make_unique<std::vector<extension::DriftChamberDigi>>()}; ///< Relation buffer for read / write
  podio::UVecPtr<edm4hep::SimTrackerHit> m_rel_sim{std::make_unique<std::vector<edm4hep::SimTrackerHit>>()}; ///< Relation buffer for read / write

  // members to handle vector members

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<MCRecoDriftChamberDigiAssociationDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociationDATA_H
#define EXTENSION_MCRecoDriftChamberDigiAssociationDATA_H


namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class MCRecoDriftChamberDigiAssociationData
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MCRecoDriftChamberDigiAssociationData {
public:
  float weight{}; ///< weight of this association

};


namespace v1 {
using MCRecoDriftChamberDigiAssociationData = extension::MCRecoDriftChamberDigiAssociationData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociationOBJ_H
#define EXTENSION_MCRecoDriftChamberDigiAssociationOBJ_H

// data model specific includes
#include "extension/MCRecoDriftChamberDigiAssociationData.h"

#include "podio/ObjectID.h"
#include <memory>
// forward declarations
namespace extension {
class DriftChamberDigi;
}
namespace edm4hep {
class SimTrackerHit;
}


namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class MCRecoDriftChamberDigiAssociation;

class MCRecoDriftChamberDigiAssociationObj {
public:
  /// constructor
  MCRecoDriftChamberDigiAssociationObj();
  /// copy constructor (does a deep-copy of relation containers)
  MCRecoDriftChamberDigiAssociationObj(const MCRecoDriftChamberDigiAssociationObj&);
  /// constructor from ObjectID and MCRecoDriftChamberDigiAssociationData
  /// does not initialize the internal relation containers
  MCRecoDriftChamberDigiAssociationObj(const podio::ObjectID& id, const MCRecoDriftChamberDigiAssociationData& data);
  /// No assignment operator
  MCRecoDriftChamberDigiAssociationObj& operator=(const MCRecoDriftChamberDigiAssociationObj&) = delete;
  virtual ~MCRecoDriftChamberDigiAssociationObj();

public:
  podio::ObjectID id{};
  MCRecoDriftChamberDigiAssociationData data;
  std::unique_ptr<extension::DriftChamberDigi> m_digi{nullptr};
  std::unique_ptr<edm4hep::SimTrackerHit> m_sim{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiV2Association_H
#define EXTENSION_MCRecoDriftChamberDigiV2Association_H

#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"


#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class MCRecoDriftChamberDigiV2AssociationCollection;
class DriftChamberDigiV2;
class MutableDriftChamberDigiV2;
}
namespace edm4hep {
class SimTrackerHit;
class MutableSimTrackerHit;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::MCRecoDriftChamberDigiV2Association& obj);
};

namespace extension {

class MutableMCRecoDriftChamberDigiV2Association;
class MCRecoDriftChamberDigiV2AssociationCollection;
class MCRecoDriftChamberDigiV2AssociationCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class MCRecoDriftChamberDigiV2Association
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MCRecoDriftChamberDigiV2Association {

  friend class MutableMCRecoDriftChamberDigiV2Association;
  friend class MCRecoDriftChamberDigiV2AssociationCollection;
  friend class extension::MCRecoDriftChamberDigiV2AssociationCollectionData;
  friend class MCRecoDriftChamberDigiV2AssociationCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const MCRecoDriftChamberDigiV2Association & obj);

public:
  using mutable_type = MutableMCRecoDriftChamberDigiV2Association;
  using collection_type = MCRecoDriftChamberDigiV2AssociationCollection;

  /// default constructor
  MCRecoDriftChamberDigiV2Association() = default;

  /// Constructor initializing all members
  MCRecoDriftChamberDigiV2Association(const float weight);

  /// copy constructor
  MCRecoDriftChamberDigiV2Association(const MCRecoDriftChamberDigiV2Association& other) noexcept = default;

  /// move constructor
  MCRecoDriftChamberDigiV2Association(MCRecoDriftChamberDigiV2Association&& other) noexcept = default;

  /// copy-assignment operator
  MCRecoDriftChamberDigiV2Association& operator=(MCRecoDriftChamberDigiV2Association other) &; // Rebind this to other's internal object
  MCRecoDriftChamberDigiV2Association& operator=(MCRecoDriftChamberDigiV2Association other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableMCRecoDriftChamberDigiV2Association clone(bool cloneRelations=true) const;

  /// destructor
  ~MCRecoDriftChamberDigiV2Association() = default;

  /// converting constructor from mutable object
  MCRecoDriftChamberDigiV2Association(const MutableMCRecoDriftChamberDigiV2Association& other);

  static MCRecoDriftChamberDigiV2Association makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::MCRecoDriftChamberDigiV2Association";

  /// Access the weight of this association
  float getWeight() const;


  /// Access the reference to the digitized hit
  const extension::DriftChamberDigiV2 getDigi() const;
  /// Access the reference to the simulated hit
  const edm4hep::SimTrackerHit getSim() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from MCRecoDriftChamberDigiV2AssociationObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr}; }

  constexpr bool operator==(const MCRecoDriftChamberDigiV2Association& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableMCRecoDriftChamberDigiV2Association& other) const;

  constexpr bool operator!=(const MCRecoDriftChamberDigiV2Association& other) const { return !(*this == other); }
  bool operator!=(const MutableMCRecoDriftChamberDigiV2Association& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MCRecoDriftChamberDigiV2Association& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MCRecoDriftChamberDigiV2Association>;

  friend void swap(MCRecoDriftChamberDigiV2Association& a, MCRecoDriftChamberDigiV2Association& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing MCRecoDriftChamberDigiV2AssociationObj
  explicit MCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> obj);
  MCRecoDriftChamberDigiV2Association(MCRecoDriftChamberDigiV2AssociationObj* obj);

  podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> m_obj{new MCRecoDriftChamberDigiV2AssociationObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiV2Association& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiV2Association& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MCRecoDriftChamberDigiV2Association> {
  std::size_t operator()(const extension::MCRecoDriftChamberDigiV2Association& obj) const {
    return std::hash<extension::MCRecoDriftChamberDigiV2AssociationObj*>{}(obj.m_obj.get());
  }
};


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiV2Association::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiV2Association::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiV2Association_CollectionData_H
#define EXTENSION_MCRecoDriftChamberDigiV2Association_CollectionData_H

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiV2AssociationData.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"
#include "extension/MCRecoDriftChamberDigiV2Association.h"
#include "extension/DriftChamberDigiV2.h"
#include "edm4hep/SimTrackerHit.h"

// schema evolution specific includes

// podio specific includes
#include "podio/CollectionBuffers.h"
#include "podio/ICollectionProvider.h"

#include <deque>
#include <memory>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


using MCRecoDriftChamberDigiV2AssociationObjPointerContainer = std::deque<MCRecoDriftChamberDigiV2AssociationObj*>;
using MCRecoDriftChamberDigiV2AssociationDataContainer = std::vector<MCRecoDriftChamberDigiV2AssociationData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class MCRecoDriftChamberDigiV2AssociationCollectionData {
public:
  /**
   * The Objs of this collection
   */
  MCRecoDriftChamberDigiV2AssociationObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  MCRecoDriftChamberDigiV2AssociationCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  MCRecoDriftChamberDigiV2AssociationCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  MCRecoDriftChamberDigiV2AssociationCollectionData(const MCRecoDriftChamberDigiV2AssociationCollectionData&) = delete;
  MCRecoDriftChamberDigiV2AssociationCollectionData& operator=(const MCRecoDriftChamberDigiV2AssociationCollectionData&) = delete;
  MCRecoDriftChamberDigiV2AssociationCollectionData(MCRecoDriftChamberDigiV2AssociationCollectionData&& other) = default;
  MCRecoDriftChamberDigiV2AssociationCollectionData& operator=(MCRecoDriftChamberDigiV2AssociationCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~MCRecoDriftChamberDigiV2AssociationCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();


  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations
  podio::UVecPtr<extension::DriftChamberDigiV2> m_rel_digi{std::make_unique<std::vector<extension::DriftChamberDigiV2>>()}; ///< Relation buffer for read / write
  podio::UVecPtr<edm4hep::SimTrackerHit> m_rel_sim{std::make_unique<std::vector<edm4hep::SimTrackerHit>>()}; ///< Relation buffer for read / write

  // members to handle vector members

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<MCRecoDriftChamberDigiV2AssociationDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiV2AssociationDATA_H
#define EXTENSION_MCRecoDriftChamberDigiV2AssociationDATA_H


namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class MCRecoDriftChamberDigiV2AssociationData
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MCRecoDriftChamberDigiV2AssociationData {
public:
  float weight{}; ///< weight of this association

};


namespace v1 {
using MCRecoDriftChamberDigiV2AssociationData = extension::MCRecoDriftChamberDigiV2AssociationData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiV2AssociationOBJ_H
#define EXTENSION_MCRecoDriftChamberDigiV2AssociationOBJ_H

// data model specific includes
#include "extension/MCRecoDriftChamberDigiV2AssociationData.h"

#include "podio/ObjectID.h"
#include <memory>
// forward declarations
namespace extension {
class DriftChamberDigiV2;
}
namespace edm4hep {
class SimTrackerHit;
}


namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class MCRecoDriftChamberDigiV2Association;

class MCRecoDriftChamberDigiV2AssociationObj {
public:
  /// constructor
  MCRecoDriftChamberDigiV2AssociationObj();
  /// copy constructor (does a deep-copy of relation containers)
  MCRecoDriftChamberDigiV2AssociationObj(const MCRecoDriftChamberDigiV2AssociationObj&);
  /// constructor from ObjectID and MCRecoDriftChamberDigiV2AssociationData
  /// does not initialize the internal relation containers
  MCRecoDriftChamberDigiV2AssociationObj(const podio::ObjectID& id, const MCRecoDriftChamberDigiV2AssociationData& data);
  /// No assignment operator
  MCRecoDriftChamberDigiV2AssociationObj& operator=(const MCRecoDriftChamberDigiV2AssociationObj&) = delete;
  virtual ~MCRecoDriftChamberDigiV2AssociationObj();

public:
  podio::ObjectID id{};
  MCRecoDriftChamberDigiV2AssociationData data;
  std::unique_ptr<extension::DriftChamberDigiV2> m_digi{nullptr};
  std::unique_ptr<edm4hep::SimTrackerHit> m_sim{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableDriftChamberDigi_H
#define EXTENSION_MutableDriftChamberDigi_H

#include "extension/DriftChamberDigiObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/DriftChamberDigi.h"

#include "edm4hep/Vector3d.h"
#include <cstdint>

#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class DriftChamberDigiCollection;
class TrackerHit;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableDriftChamberDigi
 *  Drift chamber digitized hit (before tracking) in global coordinates. Assumes that the hits are radially in the middle of the cells
 *  @author: B. Francois, CERN
 */
class MutableDriftChamberDigi {

  friend class DriftChamberDigiCollection;
  friend class DriftChamberDigiMutableCollectionIterator;
  friend class DriftChamberDigi;

public:
  using object_type = DriftChamberDigi;
  using collection_type = DriftChamberDigiCollection;

  /// default constructor
  MutableDriftChamberDigi() = default;

  /// Constructor initializing all members
  MutableDriftChamberDigi(const std::uint64_t cellID, const edm4hep::Vector3d& leftPosition, const edm4hep::Vector3d& rightPosition, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount, const std::int32_t type, const std::int32_t quality, const edm4hep::Vector3d& position);

  /// copy constructor
  MutableDriftChamberDigi(const MutableDriftChamberDigi& other) noexcept = default;

  /// move constructor
  MutableDriftChamberDigi(MutableDriftChamberDigi&& other) noexcept = default;

  /// copy-assignment operator
  MutableDriftChamberDigi& operator=(MutableDriftChamberDigi other) &; // Rebind this to other's internal object
  MutableDriftChamberDigi& operator=(MutableDriftChamberDigi other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigi clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableDriftChamberDigi() = default;


public:

  /// Access the ID of the wire that created this hit
  std::uint64_t getCellID() const;

  /// Access the position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]
  const edm4hep::Vector3d& getLeftPosition() const;

  /// Access the position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]
  const edm4hep::Vector3d& getRightPosition() const;

  /// Access the time of the hit [ns].
  float getTime() const;

  /// Access the energy deposited on the hit [GeV].
  float getEDep() const;

  /// Access the error measured on eDep [GeV].
  float getEDepError() const;

  /// Access the number of clusters associated to this hit
  std::uint32_t getClusterCount() const;

  /// Access the type of the raw data hit
  std::int32_t getType() const;

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const;

  /// Access the hit position [mm]
  const edm4hep::Vector3d& getPosition() const;



  /// Set the ID of the wire that created this hit
  void setCellID(const std::uint64_t cellID);
  /// Get mutable reference to ID of the wire that created this hit
  std::uint64_t& getCellID();

  /// Set the position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]
  void setLeftPosition(const edm4hep::Vector3d& leftPosition);
  /// Get mutable reference to position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]
  edm4hep::Vector3d& getLeftPosition();

  /// Set the position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]
  void setRightPosition(const edm4hep::Vector3d& rightPosition);
  /// Get mutable reference to position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]
  edm4hep::Vector3d& getRightPosition();

  /// Set the time of the hit [ns].
  void setTime(const float time);
  /// Get mutable reference to time of the hit [ns].
  float& getTime();

  /// Set the energy deposited on the hit [GeV].
  void setEDep(const float eDep);
  /// Get mutable reference to energy deposited on the hit [GeV].
  float& getEDep();

  /// Set the error measured on eDep [GeV].
  void setEDepError(const float eDepError);
  /// Get mutable reference to error measured on eDep [GeV].
  float& getEDepError();

  /// Set the number of clusters associated to this hit
  void setClusterCount(const std::uint32_t clusterCount);
  /// Get mutable reference to number of clusters associated to this hit
  std::uint32_t& getClusterCount();

  /// Set the type of the raw data hit
  void setType(const std::int32_t type);
  /// Get mutable reference to type of the raw data hit
  std::int32_t& getType();

  /// Set the quality bit flag of the hit
  void setQuality(const std::int32_t quality);
  /// Get mutable reference to quality bit flag of the hit
  std::int32_t& getQuality();

  /// Set the hit position [mm]
  void setPosition(const edm4hep::Vector3d& position);
  /// Get mutable reference to hit position [mm]
  edm4hep::Vector3d& getPosition();






  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiObj>{nullptr}; }

  constexpr bool operator==(const MutableDriftChamberDigi& other) const { return m_obj == other.m_obj; }
  bool operator==(const DriftChamberDigi& other) const;

  constexpr bool operator!=(const MutableDriftChamberDigi& other) const { return !(*this == other); }
  bool operator!=(const DriftChamberDigi& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableDriftChamberDigi& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableDriftChamberDigi>;

  friend void swap(MutableDriftChamberDigi& a, MutableDriftChamberDigi& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiObj
  explicit MutableDriftChamberDigi(podio::utils::MaybeSharedPtr<DriftChamberDigiObj> obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiObj> m_obj{new DriftChamberDigiObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableDriftChamberDigi& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableDriftChamberDigi> {
  std::size_t operator()(const extension::MutableDriftChamberDigi& obj) const {
    return std::hash<extension::DriftChamberDigiObj*>{}(obj.m_obj.get());
  }
};


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableDriftChamberDigiLocal_H
#define EXTENSION_MutableDriftChamberDigiLocal_H

#include "extension/DriftChamberDigiLocalObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/DriftChamberDigiLocal.h"

#include <cstdint>

#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class DriftChamberDigiLocalCollection;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableDriftChamberDigiLocal
 *  Drift chamber digitized hit (before tracking) in local coordinates
 *  @author: B. Francois, CERN
 */
class MutableDriftChamberDigiLocal {

  friend class DriftChamberDigiLocalCollection;
  friend class DriftChamberDigiLocalMutableCollectionIterator;
  friend class DriftChamberDigiLocal;

public:
  using object_type = DriftChamberDigiLocal;
  using collection_type = DriftChamberDigiLocalCollection;

  /// default constructor
  MutableDriftChamberDigiLocal() = default;

  /// Constructor initializing all members
  MutableDriftChamberDigiLocal(const std::uint64_t cellID, const float distanceToWire, const float zPositionAlongWire, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount);

  /// copy constructor
  MutableDriftChamberDigiLocal(const MutableDriftChamberDigiLocal& other) noexcept = default;

  /// move constructor
  MutableDriftChamberDigiLocal(MutableDriftChamberDigiLocal&& other) noexcept = default;

  /// copy-assignment operator
  MutableDriftChamberDigiLocal& operator=(MutableDriftChamberDigiLocal other) &; // Rebind this to other's internal object
  MutableDriftChamberDigiLocal& operator=(MutableDriftChamberDigiLocal other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigiLocal clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableDriftChamberDigiLocal() = default;


public:

  /// Access the ID of the wire that created this hit
  std::uint64_t getCellID() const;

  /// Access the smeared distance of closest approach between the wire and the hit [mm]
  float getDistanceToWire() const;

  /// Access the smeared z position in the local wire coordinate system [mm]
  float getZPositionAlongWire() const;

  /// Access the time of the hit [ns].
  float getTime() const;

  /// Access the energy deposited on the hit [GeV].
  float getEDep() const;

  /// Access the error measured on eDep [GeV].
  float getEDepError() const;

  /// Access the number of clusters associated to this hit
  std::uint32_t getClusterCount() const;



  /// Set the ID of the wire that created this hit
  void setCellID(const std::uint64_t cellID);
  /// Get mutable reference to ID of the wire that created this hit
  std::uint64_t& getCellID();

  /// Set the smeared distance of closest approach between the wire and the hit [mm]
  void setDistanceToWire(const float distanceToWire);
  /// Get mutable reference to smeared distance of closest approach between the wire and the hit [mm]
  float& getDistanceToWire();

  /// Set the smeared z position in the local wire coordinate system [mm]
  void setZPositionAlongWire(const float zPositionAlongWire);
  /// Get mutable reference to smeared z position in the local wire coordinate system [mm]
  float& getZPositionAlongWire();

  /// Set the time of the hit [ns].
  void setTime(const float time);
  /// Get mutable reference to time of the hit [ns].
  float& getTime();

  /// Set the energy deposited on the hit [GeV].
  void setEDep(const float eDep);
  /// Get mutable reference to energy deposited on the hit [GeV].
  float& getEDep();

  /// Set the error measured on eDep [GeV].
  void setEDepError(const float eDepError);
  /// Get mutable reference to error measured on eDep [GeV].
  float& getEDepError();

  /// Set the number of clusters associated to this hit
  void setClusterCount(const std::uint32_t clusterCount);
  /// Get mutable reference to number of clusters associated to this hit
  std::uint32_t& getClusterCount();






  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiLocalObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>{nullptr}; }

  constexpr bool operator==(const MutableDriftChamberDigiLocal& other) const { return m_obj == other.m_obj; }
  bool operator==(const DriftChamberDigiLocal& other) const;

  constexpr bool operator!=(const MutableDriftChamberDigiLocal& other) const { return !(*this == other); }
  bool operator!=(const DriftChamberDigiLocal& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableDriftChamberDigiLocal& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableDriftChamberDigiLocal>;

  friend void swap(MutableDriftChamberDigiLocal& a, MutableDriftChamberDigiLocal& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiLocalObj
  explicit MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> m_obj{new DriftChamberDigiLocalObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableDriftChamberDigiLocal& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableDriftChamberDigiLocal> {
  std::size_t operator()(const extension::MutableDriftChamberDigiLocal& obj) const {
    return std::hash<extension::DriftChamberDigiLocalObj*>{}(obj.m_obj.get());
  }
};


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableDriftChamberDigiV2_H
#define EXTENSION_MutableDriftChamberDigiV2_H

#include "extension/DriftChamberDigiV2Obj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/DriftChamberDigiV2.h"

#include "edm4hep/Vector3d.h"
#include "podio/RelationRange.h"
#include <cstdint>
#include <vector>

#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class DriftChamberDigiV2Collection;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableDriftChamberDigiV2
 *  Digitized hit (before tracking) for Drift Chamber v2 (requires data extension).
 *  @author: A. Tolosa-Delgado, B. Francois, CERN
 */
class MutableDriftChamberDigiV2 {

  friend class DriftChamberDigiV2Collection;
  friend class DriftChamberDigiV2MutableCollectionIterator;
  friend class DriftChamberDigiV2;

public:
  using object_type = DriftChamberDigiV2;
  using collection_type = DriftChamberDigiV2Collection;

  /// default constructor
  MutableDriftChamberDigiV2() = default;

  /// Constructor initializing all members
  MutableDriftChamberDigiV2(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const edm4hep::Vector3d& position, const edm4hep::Vector3d& directionSW, const float distanceToWire, const std::uint32_t nCluster);

  /// copy constructor
  MutableDriftChamberDigiV2(const MutableDriftChamberDigiV2& other) noexcept = default;

  /// move constructor
  MutableDriftChamberDigiV2(MutableDriftChamberDigiV2&& other) noexcept = default;

  /// copy-assignment operator
  MutableDriftChamberDigiV2& operator=(MutableDriftChamberDigiV2 other) &; // Rebind this to other's internal object
  MutableDriftChamberDigiV2& operator=(MutableDriftChamberDigiV2 other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigiV2 clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableDriftChamberDigiV2() = default;


public:

  /// Access the ID of the sensor that created this hit
  std::uint64_t getCellID() const;

  /// Access the type of the raw data hit
  std::int32_t getType() const;

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const;

  /// Access the time of the hit [ns]
  float getTime() const;

  /// Access the energy deposited on the hit [GeV]
  float getEDep() const;

  /// Access the error measured on eDep [GeV]
  float getEDepError() const;

  /// Access the point on the sensitive wire (SW) which is closest to the simhit [mm]
  const edm4hep::Vector3d& getPosition() const;

  /// Access the direction of SW
  const edm4hep::Vector3d& getDirectionSW() const;

  /// Access the distance hit-wire [mm]
  float getDistanceToWire() const;

  /// Access the number of clusters associated to this hit
  std::uint32_t getNCluster() const;



  /// Set the ID of the sensor that created this hit
  void setCellID(const std::uint64_t cellID);
  /// Get mutable reference to ID of the sensor that created this hit
  std::uint64_t& getCellID();

  /// Set the type of the raw data hit
  void setType(const std::int32_t type);
  /// Get mutable reference to type of the raw data hit
  std::int32_t& getType();

  /// Set the quality bit flag of the hit
  void setQuality(const std::int32_t quality);
  /// Get mutable reference to quality bit flag of the hit
  std::int32_t& getQuality();

  /// Set the time of the hit [ns]
  void setTime(const float time);
  /// Get mutable reference to time of the hit [ns]
  float& getTime();

  /// Set the energy deposited on the hit [GeV]
  void setEDep(const float eDep);
  /// Get mutable reference to energy deposited on the hit [GeV]
  float& getEDep();

  /// Set the error measured on eDep [GeV]
  void setEDepError(const float eDepError);
  /// Get mutable reference to error measured on eDep [GeV]
  float& getEDepError();

  /// Set the point on the sensitive wire (SW) which is closest to the simhit [mm]
  void setPosition(const edm4hep::Vector3d& position);
  /// Get mutable reference to point on the sensitive wire (SW) which is closest to the simhit [mm]
  edm4hep::Vector3d& getPosition();

  /// Set the direction of SW
  void setDirectionSW(const edm4hep::Vector3d& directionSW);
  /// Get mutable reference to direction of SW
  edm4hep::Vector3d& getDirectionSW();

  /// Set the distance hit-wire [mm]
  void setDistanceToWire(const float distanceToWire);
  /// Get mutable reference to distance hit-wire [mm]
  float& getDistanceToWire();

  /// Set the number of clusters associated to this hit
  void setNCluster(const std::uint32_t nCluster);
  /// Get mutable reference to number of clusters associated to this hit
  std::uint32_t& getNCluster();



  void addToNElectrons(const std::uint16_t&);
  std::size_t nElectrons_size() const;
  std::uint16_t getNElectrons(std::size_t) const;
  std::vector<std::uint16_t>::const_iterator nElectrons_begin() const;
  std::vector<std::uint16_t>::const_iterator nElectrons_end() const;
  podio::RelationRange<std::uint16_t> getNElectrons() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiV2Obj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr}; }

  constexpr bool operator==(const MutableDriftChamberDigiV2& other) const { return m_obj == other.m_obj; }
  bool operator==(const DriftChamberDigiV2& other) const;

  constexpr bool operator!=(const MutableDriftChamberDigiV2& other) const { return !(*this == other); }
  bool operator!=(const DriftChamberDigiV2& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableDriftChamberDigiV2& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableDriftChamberDigiV2>;

  friend void swap(MutableDriftChamberDigiV2& a, MutableDriftChamberDigiV2& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiV2Obj
  explicit MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj> obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj> m_obj{new DriftChamberDigiV2Obj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableDriftChamberDigiV2& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableDriftChamberDigiV2> {
  std::size_t operator()(const extension::MutableDriftChamberDigiV2& obj) const {
    return std::hash<extension::DriftChamberDigiV2Obj*>{}(obj.m_obj.get());
  }
};


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableMCRecoDriftChamberDigiAssociation_H
#define EXTENSION_MutableMCRecoDriftChamberDigiAssociation_H

#include "extension/MCRecoDriftChamberDigiAssociationObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/MCRecoDriftChamberDigiAssociation.h"


#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class MCRecoDriftChamberDigiAssociationCollection;
class DriftChamberDigi;
class MutableDriftChamberDigi;
}
namespace edm4hep {
class SimTrackerHit;
class MutableSimTrackerHit;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableMCRecoDriftChamberDigiAssociation
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MutableMCRecoDriftChamberDigiAssociation {

  friend class MCRecoDriftChamberDigiAssociationCollection;
  friend class MCRecoDriftChamberDigiAssociationMutableCollectionIterator;
  friend class MCRecoDriftChamberDigiAssociation;

public:
  using object_type = MCRecoDriftChamberDigiAssociation;
  using collection_type = MCRecoDriftChamberDigiAssociationCollection;

  /// default constructor
  MutableMCRecoDriftChamberDigiAssociation() = default;

  /// Constructor initializing all members
  MutableMCRecoDriftChamberDigiAssociation(const float weight);

  /// copy constructor
  MutableMCRecoDriftChamberDigiAssociation(const MutableMCRecoDriftChamberDigiAssociation& other) noexcept = default;

  /// move constructor
  MutableMCRecoDriftChamberDigiAssociation(MutableMCRecoDriftChamberDigiAssociation&& other) noexcept = default;

  /// copy-assignment operator
  MutableMCRecoDriftChamberDigiAssociation& operator=(MutableMCRecoDriftChamberDigiAssociation other) &; // Rebind this to other's internal object
  MutableMCRecoDriftChamberDigiAssociation& operator=(MutableMCRecoDriftChamberDigiAssociation other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableMCRecoDriftChamberDigiAssociation clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableMCRecoDriftChamberDigiAssociation() = default;


public:

  /// Access the weight of this association
  float getWeight() const;


  /// Access the reference to the digitized hit
  const extension::DriftChamberDigi getDigi() const;
  /// Access the reference to the simulated hit
  const edm4hep::SimTrackerHit getSim() const;

  /// Set the weight of this association
  void setWeight(const float weight);
  /// Get mutable reference to weight of this association
  float& getWeight();


  /// Set the reference to the digitized hit
  void setDigi(const extension::DriftChamberDigi& value);
  /// Set the reference to the simulated hit
  void setSim(const edm4hep::SimTrackerHit& value);




  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from MCRecoDriftChamberDigiAssociationObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr}; }

  constexpr bool operator==(const MutableMCRecoDriftChamberDigiAssociation& other) const { return m_obj == other.m_obj; }
  bool operator==(const MCRecoDriftChamberDigiAssociation& other) const;

  constexpr bool operator!=(const MutableMCRecoDriftChamberDigiAssociation& other) const { return !(*this == other); }
  bool operator!=(const MCRecoDriftChamberDigiAssociation& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableMCRecoDriftChamberDigiAssociation& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableMCRecoDriftChamberDigiAssociation>;

  friend void swap(MutableMCRecoDriftChamberDigiAssociation& a, MutableMCRecoDriftChamberDigiAssociation& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing MCRecoDriftChamberDigiAssociationObj
  explicit MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> obj);

  podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> m_obj{new MCRecoDriftChamberDigiAssociationObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableMCRecoDriftChamberDigiAssociation& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableMCRecoDriftChamberDigiAssociation> {
  std::size_t operator()(const extension::MutableMCRecoDriftChamberDigiAssociation& obj) const {
    return std::hash<extension::MCRecoDriftChamberDigiAssociationObj*>{}(obj.m_obj.get());
  }
};


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableMCRecoDriftChamberDigiV2Association_H
#define EXTENSION_MutableMCRecoDriftChamberDigiV2Association_H

#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/MCRecoDriftChamberDigiV2Association.h"


#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class MCRecoDriftChamberDigiV2AssociationCollection;
class DriftChamberDigiV2;
class MutableDriftChamberDigiV2;
}
namespace edm4hep {
class SimTrackerHit;
class MutableSimTrackerHit;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableMCRecoDriftChamberDigiV2Association
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MutableMCRecoDriftChamberDigiV2Association {

  friend class MCRecoDriftChamberDigiV2AssociationCollection;
  friend class MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator;
  friend class MCRecoDriftChamberDigiV2Association;

public:
  using object_type = MCRecoDriftChamberDigiV2Association;
  using collection_type = MCRecoDriftChamberDigiV2AssociationCollection;

  /// default constructor
  MutableMCRecoDriftChamberDigiV2Association() = default;

  /// Constructor initializing all members
  MutableMCRecoDriftChamberDigiV2Association(const float weight);

  /// copy constructor
  MutableMCRecoDriftChamberDigiV2Association(const MutableMCRecoDriftChamberDigiV2Association& other) noexcept = default;

  /// move constructor
  MutableMCRecoDriftChamberDigiV2Association(MutableMCRecoDriftChamberDigiV2Association&& other) noexcept = default;

  /// copy-assignment operator
  MutableMCRecoDriftChamberDigiV2Association& operator=(MutableMCRecoDriftChamberDigiV2Association other) &; // Rebind this to other's internal object
  MutableMCRecoDriftChamberDigiV2Association& operator=(MutableMCRecoDriftChamberDigiV2Association other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableMCRecoDriftChamberDigiV2Association clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableMCRecoDriftChamberDigiV2Association() = default;


public:

  /// Access the weight of this association
  float getWeight() const;


  /// Access the reference to the digitized hit
  const extension::DriftChamberDigiV2 getDigi() const;
  /// Access the reference to the simulated hit
  const edm4hep::SimTrackerHit getSim() const;

  /// Set the weight of this association
  void setWeight(const float weight);
  /// Get mutable reference to weight of this association
  float& getWeight();


  /// Set the reference to the digitized hit
  void setDigi(const extension::DriftChamberDigiV2& value);
  /// Set the reference to the simulated hit
  void setSim(const edm4hep::SimTrackerHit& value);




  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from MCRecoDriftChamberDigiV2AssociationObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr}; }

  constexpr bool operator==(const MutableMCRecoDriftChamberDigiV2Association& other) const { return m_obj == other.m_obj; }
  bool operator==(const MCRecoDriftChamberDigiV2Association& other) const;

  constexpr bool operator!=(const MutableMCRecoDriftChamberDigiV2Association& other) const { return !(*this == other); }
  bool operator!=(const MCRecoDriftChamberDigiV2Association& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableMCRecoDriftChamberDigiV2Association& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableMCRecoDriftChamberDigiV2Association>;

  friend void swap(MutableMCRecoDriftChamberDigiV2Association& a, MutableMCRecoDriftChamberDigiV2Association& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing MCRecoDriftChamberDigiV2AssociationObj
  explicit MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> obj);

  podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> m_obj{new MCRecoDriftChamberDigiV2AssociationObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableMCRecoDriftChamberDigiV2Association& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableMCRecoDriftChamberDigiV2Association> {
  std::size_t operator()(const extension::MutableMCRecoDriftChamberDigiV2Association& obj) const {
    return std::hash<extension::MCRecoDriftChamberDigiV2AssociationObj*>{}(obj.m_obj.get());
  }
};


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableSenseWireHit_H
#define EXTENSION_MutableSenseWireHit_H

#include "extension/SenseWireHitObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/SenseWireHit.h"

#include "edm4hep/Vector3d.h"
#include "podio/RelationRange.h"
#include <cstdint>
#include <vector>

#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class SenseWireHitCollection;
class TrackerHit;
class TrackerHit_dev;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableSenseWireHit
 *  Sense wire hit, before the hit position degeneracy is resolved. The circle representing possible positions is parametrized with its center, radius and normal vector (given by the wire direction).
 *  @author: EDM4hep authors
 */
class MutableSenseWireHit {

  friend class SenseWireHitCollection;
  friend class SenseWireHitMutableCollectionIterator;
  friend class SenseWireHit;

public:
  using object_type = SenseWireHit;
  using collection_type = SenseWireHitCollection;

  /// default constructor
  MutableSenseWireHit() = default;

  /// Constructor initializing all members
  MutableSenseWireHit(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const float wireStereoAngle, const float wireAzimuthalAngle, const edm4hep::Vector3d& position, const double positionAlongWireError, const float distanceToWire, const float distanceToWireError);

  /// copy constructor
  MutableSenseWireHit(const MutableSenseWireHit& other) noexcept = default;

  /// move constructor
  MutableSenseWireHit(MutableSenseWireHit&& other) noexcept = default;

  /// copy-assignment operator
  MutableSenseWireHit& operator=(MutableSenseWireHit other) &; // Rebind this to other's internal object
  MutableSenseWireHit& operator=(MutableSenseWireHit other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableSenseWireHit clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableSenseWireHit() = default;


public:

  /// Access the ID of the sensor that created this hit
  std::uint64_t getCellID() const;

  /// Access the type of the raw data hit
  std::int32_t getType() const;

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const;

  /// Access the time of the hit [ns]
  float getTime() const;

  /// Access the energy deposited by the hit [GeV]
  float getEDep() const;

  /// Access the error on eDep [GeV]
  float getEDepError() const;

  /// Access the angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
  float getWireStereoAngle() const;

  /// Access the azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
  float getWireAzimuthalAngle() const;

  /// Access the point on the sense wire which is closest to the hit (center of the circle) [mm]
  const edm4hep::Vector3d& getPosition() const;

  /// Access the error on the hit position along the wire direction [mm]
  double getPositionAlongWireError() const;

  /// Access the distance between the hit and the wire (radius of the circle) [mm]
  float getDistanceToWire() const;

  /// Access the error on distanceToWire [mm]
  float getDistanceToWireError() const;



  /// Set the ID of the sensor that created this hit
  void setCellID(const std::uint64_t cellID);
  /// Get mutable reference to ID of the sensor that created this hit
  std::uint64_t& getCellID();

  /// Set the type of the raw data hit
  void setType(const std::int32_t type);
  /// Get mutable reference to type of the raw data hit
  std::int32_t& getType();

  /// Set the quality bit flag of the hit
  void setQuality(const std::int32_t quality);
  /// Get mutable reference to quality bit flag of the hit
  std::int32_t& getQuality();

  /// Set the time of the hit [ns]
  void setTime(const float time);
  /// Get mutable reference to time of the hit [ns]
  float& getTime();

  /// Set the energy deposited by the hit [GeV]
  void setEDep(const float eDep);
  /// Get mutable reference to energy deposited by the hit [GeV]
  float& getEDep();

  /// Set the error on eDep [GeV]
  void setEDepError(const float eDepError);
  /// Get mutable reference to error on eDep [GeV]
  float& getEDepError();

  /// Set the angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
  void setWireStereoAngle(const float wireStereoAngle);
  /// Get mutable reference to angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
  float& getWireStereoAngle();

  /// Set the azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
  void setWireAzimuthalAngle(const float wireAzimuthalAngle);
  /// Get mutable reference to azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
  float& getWireAzimuthalAngle();

  /// Set the point on the sense wire which is closest to the hit (center of the circle) [mm]
  void setPosition(const edm4hep::Vector3d& position);
  /// Get mutable reference to point on the sense wire which is closest to the hit (center of the circle) [mm]
  edm4hep::Vector3d& getPosition();

  /// Set the error on the hit position along the wire direction [mm]
  void setPositionAlongWireError(const double positionAlongWireError);
  /// Get mutable reference to error on the hit position along the wire direction [mm]
  double& getPositionAlongWireError();

  /// Set the distance between the hit and the wire (radius of the circle) [mm]
  void setDistanceToWire(const float distanceToWire);
  /// Get mutable reference to distance between the hit and the wire (radius of the circle) [mm]
  float& getDistanceToWire();

  /// Set the error on distanceToWire [mm]
  void setDistanceToWireError(const float distanceToWireError);
  /// Get mutable reference to error on distanceToWire [mm]
  float& getDistanceToWireError();



  void addToNElectrons(const std::uint16_t&);
  std::size_t nElectrons_size() const;
  std::uint16_t getNElectrons(std::size_t) const;
  std::vector<std::uint16_t>::const_iterator nElectrons_begin() const;
  std::vector<std::uint16_t>::const_iterator nElectrons_end() const;
  podio::RelationRange<std::uint16_t> getNElectrons() const;

 /// Return the number of clusters associated to the hit
 auto getNClusters() const { return getNElectrons().size(); }
 


  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from SenseWireHitObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr}; }

  constexpr bool operator==(const MutableSenseWireHit& other) const { return m_obj == other.m_obj; }
  bool operator==(const SenseWireHit& other) const;

  constexpr bool operator!=(const MutableSenseWireHit& other) const { return !(*this == other); }
  bool operator!=(const SenseWireHit& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableSenseWireHit& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableSenseWireHit>;

  friend void swap(MutableSenseWireHit& a, MutableSenseWireHit& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing SenseWireHitObj
  explicit MutableSenseWireHit(podio::utils::MaybeSharedPtr<SenseWireHitObj> obj);

  podio::utils::MaybeSharedPtr<SenseWireHitObj> m_obj{new SenseWireHitObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableSenseWireHit& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableSenseWireHit> {
  std::size_t operator()(const extension::MutableSenseWireHit& obj) const {
    return std::hash<extension::SenseWireHitObj*>{}(obj.m_obj.get());
  }
};


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableSenseWireHitSimTrackerHitLink_H
#define EXTENSION_MutableSenseWireHitSimTrackerHitLink_H

#include "extension/SenseWireHitSimTrackerHitLinkObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/SenseWireHitSimTrackerHitLink.h"


#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class SenseWireHitSimTrackerHitLinkCollection;
class SenseWireHit;
class MutableSenseWireHit;
}
namespace edm4hep {
class SimTrackerHit;
class MutableSimTrackerHit;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableSenseWireHitSimTrackerHitLink
 *  Link between a SenseWireHit and a SimTrackerHit
 *  @author: EDM4hep authors
 */
class MutableSenseWireHitSimTrackerHitLink {

  friend class SenseWireHitSimTrackerHitLinkCollection;
  friend class SenseWireHitSimTrackerHitLinkMutableCollectionIterator;
  friend class SenseWireHitSimTrackerHitLink;

public:
  using object_type = SenseWireHitSimTrackerHitLink;
  using collection_type = SenseWireHitSimTrackerHitLinkCollection;

  /// default constructor
  MutableSenseWireHitSimTrackerHitLink() = default;

  /// Constructor initializing all members
  MutableSenseWireHitSimTrackerHitLink(const float weight);

  /// copy constructor
  MutableSenseWireHitSimTrackerHitLink(const MutableSenseWireHitSimTrackerHitLink& other) noexcept = default;

  /// move constructor
  MutableSenseWireHitSimTrackerHitLink(MutableSenseWireHitSimTrackerHitLink&& other) noexcept = default;

  /// copy-assignment operator
  MutableSenseWireHitSimTrackerHitLink& operator=(MutableSenseWireHitSimTrackerHitLink other) &; // Rebind this to other's internal object
  MutableSenseWireHitSimTrackerHitLink& operator=(MutableSenseWireHitSimTrackerHitLink other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableSenseWireHitSimTrackerHitLink clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableSenseWireHitSimTrackerHitLink() = default;


public:

  /// Access the weight of this link
  float getWeight() const;


  /// Access the reference to the SenseWireHit
  const extension::SenseWireHit getFrom() const;
  /// Access the reference to the SimTrackerHit
  const edm4hep::SimTrackerHit getTo() const;

  /// Set the weight of this link
  void setWeight(const float weight);
  /// Get mutable reference to weight of this link
  float& getWeight();


  /// Set the reference to the SenseWireHit
  void setFrom(const extension::SenseWireHit& value);
  /// Set the reference to the SimTrackerHit
  void setTo(const edm4hep::SimTrackerHit& value);




  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from SenseWireHitSimTrackerHitLinkObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr}; }

  constexpr bool operator==(const MutableSenseWireHitSimTrackerHitLink& other) const { return m_obj == other.m_obj; }
  bool operator==(const SenseWireHitSimTrackerHitLink& other) const;

  constexpr bool operator!=(const MutableSenseWireHitSimTrackerHitLink& other) const { return !(*this == other); }
  bool operator!=(const SenseWireHitSimTrackerHitLink& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableSenseWireHitSimTrackerHitLink& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableSenseWireHitSimTrackerHitLink>;

  friend void swap(MutableSenseWireHitSimTrackerHitLink& a, MutableSenseWireHitSimTrackerHitLink& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing SenseWireHitSimTrackerHitLinkObj
  explicit MutableSenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> obj);

  podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> m_obj{new SenseWireHitSimTrackerHitLinkObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableSenseWireHitSimTrackerHitLink& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableSenseWireHitSimTrackerHitLink> {
  std::size_t operator()(const extension::MutableSenseWireHitSimTrackerHitLink& obj) const {
    return std::hash<extension::SenseWireHitSimTrackerHitLinkObj*>{}(obj.m_obj.get());
  }
};


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableTrack_H
#define EXTENSION_MutableTrack_H

#include "extension/TrackObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/Track.h"

#include "extension/TrackerHit.h"
#include "podio/RelationRange.h"
#include <cstdint>
#include <vector>

#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class TrackCollection;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableTrack
 *  Reconstructed track
 *  @author: EDM4hep authors
 */
class MutableTrack {

  friend class TrackCollection;
  friend class TrackMutableCollectionIterator;
  friend class Track;

public:
  using object_type = Track;
  using collection_type = TrackCollection;

  /// default constructor
  MutableTrack() = default;

  /// Constructor initializing all members
  MutableTrack(const std::int32_t type, const float chi2, const std::int32_t ndf, const float dEdx, const float dEdxError, const float radiusOfInnermostHit);

  /// copy constructor
  MutableTrack(const MutableTrack& other) noexcept = default;

  /// move constructor
  MutableTrack(MutableTrack&& other) noexcept = default;

  /// copy-assignment operator
  MutableTrack& operator=(MutableTrack other) &; // Rebind this to other's internal object
  MutableTrack& operator=(MutableTrack other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableTrack clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableTrack() = default;


public:

  /// Access the flagword that defines the type of track.Bits 16-31 are used internally
  std::int32_t getType() const;

  /// Access the Chi^2 of the track fit
  float getChi2() const;

  /// Access the number of degrees of freedom of the track fit
  std::int32_t getNdf() const;

  /// Access the dEdx of the track
  float getDEdx() const;

  /// Access the error of dEdx
  float getDEdxError() const;

  /// Access the radius of the innermost hit that has been used in the track fit
  float getRadiusOfInnermostHit() const;



  /// Set the flagword that defines the type of track.Bits 16-31 are used internally
  void setType(const std::int32_t type);
  /// Get mutable reference to flagword that defines the type of track.Bits 16-31 are used internally
  std::int32_t& getType();

  /// Set the Chi^2 of the track fit
  void setChi2(const float chi2);
  /// Get mutable reference to Chi^2 of the track fit
  float& getChi2();

  /// Set the number of degrees of freedom of the track fit
  void setNdf(const std::int32_t ndf);
  /// Get mutable reference to number of degrees of freedom of the track fit
  std::int32_t& getNdf();

  /// Set the dEdx of the track
  void setDEdx(const float dEdx);
  /// Get mutable reference to dEdx of the track
  float& getDEdx();

  /// Set the error of dEdx
  void setDEdxError(const float dEdxError);
  /// Get mutable reference to error of dEdx
  float& getDEdxError();

  /// Set the radius of the innermost hit that has been used in the track fit
  void setRadiusOfInnermostHit(const float radiusOfInnermostHit);
  /// Get mutable reference to radius of the innermost hit that has been used in the track fit
  float& getRadiusOfInnermostHit();



  void addToTrackerHits(const extension::TrackerHit&);
  std::size_t trackerHits_size() const;
  extension::TrackerHit getTrackerHits(std::size_t) const;
  std::vector<extension::TrackerHit>::const_iterator trackerHits_begin() const;
  std::vector<extension::TrackerHit>::const_iterator trackerHits_end() const;
  podio::RelationRange<extension::TrackerHit> getTrackerHits() const;
  void addToTracks(const extension::Track&);
  std::size_t tracks_size() const;
  extension::Track getTracks(std::size_t) const;
  std::vector<extension::Track>::const_iterator tracks_begin() const;
  std::vector<extension::Track>::const_iterator tracks_end() const;
  podio::RelationRange<extension::Track> getTracks() const;
  void addToSubdetectorHitNumbers(const std::int32_t&);
  std::size_t subdetectorHitNumbers_size() const;
  std::int32_t getSubdetectorHitNumbers(std::size_t) const;
  std::vector<std::int32_t>::const_iterator subdetectorHitNumbers_begin() const;
  std::vector<std::int32_t>::const_iterator subdetectorHitNumbers_end() const;
  podio::RelationRange<std::int32_t> getSubdetectorHitNumbers() const;
  void addToTrackStates(const edm4hep::TrackState&);
  std::size_t trackStates_size() const;
  edm4hep::TrackState getTrackStates(std::size_t) const;
  std::vector<edm4hep::TrackState>::const_iterator trackStates_begin() const;
  std::vector<edm4hep::TrackState>::const_iterator trackStates_end() const;
  podio::RelationRange<edm4hep::TrackState> getTrackStates() const;
  void addToDxQuantities(const edm4hep::Quantity&);
  std::size_t dxQuantities_size() const;
  edm4hep::Quantity getDxQuantities(std::size_t) const;
  std::vector<edm4hep::Quantity>::const_iterator dxQuantities_begin() const;
  std::vector<edm4hep::Quantity>::const_iterator dxQuantities_end() const;
  podio::RelationRange<edm4hep::Quantity> getDxQuantities() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from TrackObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<TrackObj>{nullptr}; }

  constexpr bool operator==(const MutableTrack& other) const { return m_obj == other.m_obj; }
  bool operator==(const Track& other) const;

  constexpr bool operator!=(const MutableTrack& other) const { return !(*this == other); }
  bool operator!=(const Track& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableTrack& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableTrack>;

  friend void swap(MutableTrack& a, MutableTrack& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing TrackObj
  explicit MutableTrack(podio::utils::MaybeSharedPtr<TrackObj> obj);

  podio::utils::MaybeSharedPtr<TrackObj> m_obj{new TrackObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableTrack& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableTrack> {
  std::size_t operator()(const extension::MutableTrack& obj) const {
    return std::hash<extension::TrackObj*>{}(obj.m_obj.get());
  }
};


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHit_H
#define EXTENSION_SenseWireHit_H

#include "extension/SenseWireHitObj.h"

#include "edm4hep/Vector3d.h"
#include "podio/RelationRange.h"
#include <cstdint>
#include <vector>

#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class SenseWireHitCollection;
class TrackerHit;
class TrackerHit_dev;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::SenseWireHit& obj);
};

namespace extension {

class MutableSenseWireHit;
class SenseWireHitCollection;
class SenseWireHitCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class SenseWireHit
 *  Sense wire hit, before the hit position degeneracy is resolved. The circle representing possible positions is parametrized with its center, radius and normal vector (given by the wire direction).
 *  @author: EDM4hep authors
 */
class SenseWireHit {

  friend class MutableSenseWireHit;
  friend class SenseWireHitCollection;
  friend class extension::SenseWireHitCollectionData;
  friend class SenseWireHitCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const SenseWireHit & obj);

public:
  using mutable_type = MutableSenseWireHit;
  using collection_type = SenseWireHitCollection;

  /// default constructor
  SenseWireHit() = default;

  /// Constructor initializing all members
  SenseWireHit(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const float wireStereoAngle, const float wireAzimuthalAngle, const edm4hep::Vector3d& position, const double positionAlongWireError, const float distanceToWire, const float distanceToWireError);

  /// copy constructor
  SenseWireHit(const SenseWireHit& other) noexcept = default;

  /// move constructor
  SenseWireHit(SenseWireHit&& other) noexcept = default;

  /// copy-assignment operator
  SenseWireHit& operator=(SenseWireHit other) &; // Rebind this to other's internal object
  SenseWireHit& operator=(SenseWireHit other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableSenseWireHit clone(bool cloneRelations=true) const;

  /// destructor
  ~SenseWireHit() = default;

  /// converting constructor from mutable object
  SenseWireHit(const MutableSenseWireHit& other);

  static SenseWireHit makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::SenseWireHit";

  /// Access the ID of the sensor that created this hit
  std::uint64_t getCellID() const;

  /// Access the type of the raw data hit
  std::int32_t getType() const;

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const;

  /// Access the time of the hit [ns]
  float getTime() const;

  /// Access the energy deposited by the hit [GeV]
  float getEDep() const;

  /// Access the error on eDep [GeV]
  float getEDepError() const;

  /// Access the angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
  float getWireStereoAngle() const;

  /// Access the azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
  float getWireAzimuthalAngle() const;

  /// Access the point on the sense wire which is closest to the hit (center of the circle) [mm]
  const edm4hep::Vector3d& getPosition() const;

  /// Access the error on the hit position along the wire direction [mm]
  double getPositionAlongWireError() const;

  /// Access the distance between the hit and the wire (radius of the circle) [mm]
  float getDistanceToWire() const;

  /// Access the error on distanceToWire [mm]
  float getDistanceToWireError() const;



  std::size_t nElectrons_size() const;
  std::uint16_t getNElectrons(std::size_t) const;
  std::vector<std::uint16_t>::const_iterator nElectrons_begin() const;
  std::vector<std::uint16_t>::const_iterator nElectrons_end() const;
  podio::RelationRange<std::uint16_t> getNElectrons() const;

 /// Return the number of clusters associated to the hit
 auto getNClusters() const { return getNElectrons().size(); }
 

  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from SenseWireHitObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr}; }

  constexpr bool operator==(const SenseWireHit& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableSenseWireHit& other) const;

  constexpr bool operator!=(const SenseWireHit& other) const { return !(*this == other); }
  bool operator!=(const MutableSenseWireHit& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const SenseWireHit& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<SenseWireHit>;

  friend void swap(SenseWireHit& a, SenseWireHit& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing SenseWireHitObj
  explicit SenseWireHit(podio::utils::MaybeSharedPtr<SenseWireHitObj> obj);
  SenseWireHit(SenseWireHitObj* obj);

  podio::utils::MaybeSharedPtr<SenseWireHitObj> m_obj{new SenseWireHitObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const SenseWireHit& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHit& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::SenseWireHit> {
  std::size_t operator()(const extension::SenseWireHit& obj) const {
    return std::hash<extension::SenseWireHitObj*>{}(obj.m_obj.get());
  }
};


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHit::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHit::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHit_CollectionData_H
#define EXTENSION_SenseWireHit_CollectionData_H

// datamodel specific includes
#include "extension/SenseWireHitData.h"
#include "extension/SenseWireHitObj.h"
#include "extension/SenseWireHit.h"

// schema evolution specific includes

// podio specific includes
#include "podio/CollectionBuffers.h"
#include "podio/ICollectionProvider.h"

#include <deque>
#include <memory>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


using SenseWireHitObjPointerContainer = std::deque<SenseWireHitObj*>;
using SenseWireHitDataContainer = std::vector<SenseWireHitData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class SenseWireHitCollectionData {
public:
  /**
   * The Objs of this collection
   */
  SenseWireHitObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  SenseWireHitCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  SenseWireHitCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  SenseWireHitCollectionData(const SenseWireHitCollectionData&) = delete;
  SenseWireHitCollectionData& operator=(const SenseWireHitCollectionData&) = delete;
  SenseWireHitCollectionData(SenseWireHitCollectionData&& other) = default;
  SenseWireHitCollectionData& operator=(SenseWireHitCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~SenseWireHitCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();

  void createRelations(SenseWireHitObj* obj);

  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations

  // members to handle vector members
  podio::UVecPtr<std::uint16_t> m_vec_nElectrons{nullptr}; /// combined vector of all objects in collection
  std::vector<podio::UVecPtr<std::uint16_t>> m_vecs_nElectrons{}; /// pointers to individual member vectors

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<SenseWireHitDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitDATA_H
#define EXTENSION_SenseWireHitDATA_H

#include "edm4hep/Vector3d.h"
#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class SenseWireHitData
 *  Sense wire hit, before the hit position degeneracy is resolved. The circle representing possible positions is parametrized with its center, radius and normal vector (given by the wire direction).
 *  @author: EDM4hep authors
 */
class SenseWireHitData {
public:
  std::uint64_t cellID{}; ///< ID of the sensor that created this hit
  std::int32_t type{}; ///< type of the raw data hit
  std::int32_t quality{}; ///< quality bit flag of the hit
  float time{}; ///< time of the hit [ns]
  float eDep{}; ///< energy deposited by the hit [GeV]
  float eDepError{}; ///< error on eDep [GeV]
  float wireStereoAngle{}; ///< angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
  float wireAzimuthalAngle{}; ///< azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
  ::edm4hep::Vector3d position{}; ///< point on the sense wire which is closest to the hit (center of the circle) [mm]
  double positionAlongWireError{}; ///< error on the hit position along the wire direction [mm]
  float distanceToWire{}; ///< distance between the hit and the wire (radius of the circle) [mm]
  float distanceToWireError{}; ///< error on distanceToWire [mm]

  unsigned int nElectrons_begin{};
  unsigned int nElectrons_end{};
};


namespace v1 {
using SenseWireHitData = extension::SenseWireHitData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitOBJ_H
#define EXTENSION_SenseWireHitOBJ_H

// data model specific includes
#include "extension/SenseWireHitData.h"
#include <vector>

#include "podio/ObjectID.h"
#include <vector>



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class SenseWireHit;

class SenseWireHitObj {
public:
  /// constructor
  SenseWireHitObj();
  /// copy constructor (does a deep-copy of relation containers)
  SenseWireHitObj(const SenseWireHitObj&);
  /// constructor from ObjectID and SenseWireHitData
  /// does not initialize the internal relation containers
  SenseWireHitObj(const podio::ObjectID& id, const SenseWireHitData& data);
  /// No assignment operator
  SenseWireHitObj& operator=(const SenseWireHitObj&) = delete;
  virtual ~SenseWireHitObj();

public:
  podio::ObjectID id{};
  SenseWireHitData data;
  std::vector<std::uint16_t>* m_nElectrons{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLink_H
#define EXTENSION_SenseWireHitSimTrackerHitLink_H

#include "extension/SenseWireHitSimTrackerHitLinkObj.h"


#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class SenseWireHitSimTrackerHitLinkCollection;
class SenseWireHit;
class MutableSenseWireHit;
}
namespace edm4hep {
class SimTrackerHit;
class MutableSimTrackerHit;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::SenseWireHitSimTrackerHitLink& obj);
};

namespace extension {

class MutableSenseWireHitSimTrackerHitLink;
class SenseWireHitSimTrackerHitLinkCollection;
class SenseWireHitSimTrackerHitLinkCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class SenseWireHitSimTrackerHitLink
 *  Link between a SenseWireHit and a SimTrackerHit
 *  @author: EDM4hep authors
 */
class SenseWireHitSimTrackerHitLink {

  friend class MutableSenseWireHitSimTrackerHitLink;
  friend class SenseWireHitSimTrackerHitLinkCollection;
  friend class extension::SenseWireHitSimTrackerHitLinkCollectionData;
  friend class SenseWireHitSimTrackerHitLinkCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const SenseWireHitSimTrackerHitLink & obj);

public:
  using mutable_type = MutableSenseWireHitSimTrackerHitLink;
  using collection_type = SenseWireHitSimTrackerHitLinkCollection;

  /// default constructor
  SenseWireHitSimTrackerHitLink() = default;

  /// Constructor initializing all members
  SenseWireHitSimTrackerHitLink(const float weight);

  /// copy constructor
  SenseWireHitSimTrackerHitLink(const SenseWireHitSimTrackerHitLink& other) noexcept = default;

  /// move constructor
  SenseWireHitSimTrackerHitLink(SenseWireHitSimTrackerHitLink&& other) noexcept = default;

  /// copy-assignment operator
  SenseWireHitSimTrackerHitLink& operator=(SenseWireHitSimTrackerHitLink other) &; // Rebind this to other's internal object
  SenseWireHitSimTrackerHitLink& operator=(SenseWireHitSimTrackerHitLink other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableSenseWireHitSimTrackerHitLink clone(bool cloneRelations=true) const;

  /// destructor
  ~SenseWireHitSimTrackerHitLink() = default;

  /// converting constructor from mutable object
  SenseWireHitSimTrackerHitLink(const MutableSenseWireHitSimTrackerHitLink& other);

  static SenseWireHitSimTrackerHitLink makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::SenseWireHitSimTrackerHitLink";

  /// Access the weight of this link
  float getWeight() const;


  /// Access the reference to the SenseWireHit
  const extension::SenseWireHit getFrom() const;
  /// Access the reference to the SimTrackerHit
  const edm4hep::SimTrackerHit getTo() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from SenseWireHitSimTrackerHitLinkObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr}; }

  constexpr bool operator==(const SenseWireHitSimTrackerHitLink& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableSenseWireHitSimTrackerHitLink& other) const;

  constexpr bool operator!=(const SenseWireHitSimTrackerHitLink& other) const { return !(*this == other); }
  bool operator!=(const MutableSenseWireHitSimTrackerHitLink& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const SenseWireHitSimTrackerHitLink& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<SenseWireHitSimTrackerHitLink>;

  friend void swap(SenseWireHitSimTrackerHitLink& a, SenseWireHitSimTrackerHitLink& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing SenseWireHitSimTrackerHitLinkObj
  explicit SenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> obj);
  SenseWireHitSimTrackerHitLink(SenseWireHitSimTrackerHitLinkObj* obj);

  podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> m_obj{new SenseWireHitSimTrackerHitLinkObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const SenseWireHitSimTrackerHitLink& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHitSimTrackerHitLink& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::SenseWireHitSimTrackerHitLink> {
  std::size_t operator()(const extension::SenseWireHitSimTrackerHitLink& obj) const {
    return std::hash<extension::SenseWireHitSimTrackerHitLinkObj*>{}(obj.m_obj.get());
  }
};


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitSimTrackerHitLink::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitSimTrackerHitLink::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLink_CollectionData_H
#define EXTENSION_SenseWireHitSimTrackerHitLink_CollectionData_H

// datamodel specific includes
#include "extension/SenseWireHitSimTrackerHitLinkData.h"
#include "extension/SenseWireHitSimTrackerHitLinkObj.h"
#include "extension/SenseWireHitSimTrackerHitLink.h"
#include "extension/SenseWireHit.h"
#include "edm4hep/SimTrackerHit.h"

// schema evolution specific includes

// podio specific includes
#include "podio/CollectionBuffers.h"
#include "podio/ICollectionProvider.h"

#include <deque>
#include <memory>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


using SenseWireHitSimTrackerHitLinkObjPointerContainer = std::deque<SenseWireHitSimTrackerHitLinkObj*>;
using SenseWireHitSimTrackerHitLinkDataContainer = std::vector<SenseWireHitSimTrackerHitLinkData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class SenseWireHitSimTrackerHitLinkCollectionData {
public:
  /**
   * The Objs of this collection
   */
  SenseWireHitSimTrackerHitLinkObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  SenseWireHitSimTrackerHitLinkCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  SenseWireHitSimTrackerHitLinkCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  SenseWireHitSimTrackerHitLinkCollectionData(const SenseWireHitSimTrackerHitLinkCollectionData&) = delete;
  SenseWireHitSimTrackerHitLinkCollectionData& operator=(const SenseWireHitSimTrackerHitLinkCollectionData&) = delete;
  SenseWireHitSimTrackerHitLinkCollectionData(SenseWireHitSimTrackerHitLinkCollectionData&& other) = default;
  SenseWireHitSimTrackerHitLinkCollectionData& operator=(SenseWireHitSimTrackerHitLinkCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~SenseWireHitSimTrackerHitLinkCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();


  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations
  podio::UVecPtr<extension::SenseWireHit> m_rel_from{std::make_unique<std::vector<extension::SenseWireHit>>()}; ///< Relation buffer for read / write
  podio::UVecPtr<edm4hep::SimTrackerHit> m_rel_to{std::make_unique<std::vector<edm4hep::SimTrackerHit>>()}; ///< Relation buffer for read / write

  // members to handle vector members

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<SenseWireHitSimTrackerHitLinkDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLinkDATA_H
#define EXTENSION_SenseWireHitSimTrackerHitLinkDATA_H


namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class SenseWireHitSimTrackerHitLinkData
 *  Link between a SenseWireHit and a SimTrackerHit
 *  @author: EDM4hep authors
 */
class SenseWireHitSimTrackerHitLinkData {
public:
  float weight{}; ///< weight of this link

};


namespace v1 {
using SenseWireHitSimTrackerHitLinkData = extension::SenseWireHitSimTrackerHitLinkData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLinkOBJ_H
#define EXTENSION_SenseWireHitSimTrackerHitLinkOBJ_H

// data model specific includes
#include "extension/SenseWireHitSimTrackerHitLinkData.h"

#include "podio/ObjectID.h"
#include <memory>
// forward declarations
namespace extension {
class SenseWireHit;
}
namespace edm4hep {
class SimTrackerHit;
}


namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class SenseWireHitSimTrackerHitLink;

class SenseWireHitSimTrackerHitLinkObj {
public:
  /// constructor
  SenseWireHitSimTrackerHitLinkObj();
  /// copy constructor (does a deep-copy of relation containers)
  SenseWireHitSimTrackerHitLinkObj(const SenseWireHitSimTrackerHitLinkObj&);
  /// constructor from ObjectID and SenseWireHitSimTrackerHitLinkData
  /// does not initialize the internal relation containers
  SenseWireHitSimTrackerHitLinkObj(const podio::ObjectID& id, const SenseWireHitSimTrackerHitLinkData& data);
  /// No assignment operator
  SenseWireHitSimTrackerHitLinkObj& operator=(const SenseWireHitSimTrackerHitLinkObj&) = delete;
  virtual ~SenseWireHitSimTrackerHitLinkObj();

public:
  podio::ObjectID id{};
  SenseWireHitSimTrackerHitLinkData data;
  std::unique_ptr<extension::SenseWireHit> m_from{nullptr};
  std::unique_ptr<edm4hep::SimTrackerHit> m_to{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_Track_H
#define EXTENSION_Track_H

#include "extension/TrackObj.h"

#include "extension/TrackerHit.h"
#include "podio/RelationRange.h"
#include <cstdint>
#include <vector>

#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class TrackCollection;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::Track& obj);
};

namespace extension {

class MutableTrack;
class TrackCollection;
class TrackCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class Track
 *  Reconstructed track
 *  @author: EDM4hep authors
 */
class Track {

  friend class MutableTrack;
  friend class TrackCollection;
  friend class extension::TrackCollectionData;
  friend class TrackCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const Track & obj);

public:
  using mutable_type = MutableTrack;
  using collection_type = TrackCollection;

  /// default constructor
  Track() = default;

  /// Constructor initializing all members
  Track(const std::int32_t type, const float chi2, const std::int32_t ndf, const float dEdx, const float dEdxError, const float radiusOfInnermostHit);

  /// copy constructor
  Track(const Track& other) noexcept = default;

  /// move constructor
  Track(Track&& other) noexcept = default;

  /// copy-assignment operator
  Track& operator=(Track other) &; // Rebind this to other's internal object
  Track& operator=(Track other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableTrack clone(bool cloneRelations=true) const;

  /// destructor
  ~Track() = default;

  /// converting constructor from mutable object
  Track(const MutableTrack& other);

  static Track makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::Track";

  /// Access the flagword that defines the type of track.Bits 16-31 are used internally
  std::int32_t getType() const;

  /// Access the Chi^2 of the track fit
  float getChi2() const;

  /// Access the number of degrees of freedom of the track fit
  std::int32_t getNdf() const;

  /// Access the dEdx of the track
  float getDEdx() const;

  /// Access the error of dEdx
  float getDEdxError() const;

  /// Access the radius of the innermost hit that has been used in the track fit
  float getRadiusOfInnermostHit() const;



  std::size_t trackerHits_size() const;
  extension::TrackerHit getTrackerHits(std::size_t) const;
  std::vector<extension::TrackerHit>::const_iterator trackerHits_begin() const;
  std::vector<extension::TrackerHit>::const_iterator trackerHits_end() const;
  podio::RelationRange<extension::TrackerHit> getTrackerHits() const;
  std::size_t tracks_size() const;
  extension::Track getTracks(std::size_t) const;
  std::vector<extension::Track>::const_iterator tracks_begin() const;
  std::vector<extension::Track>::const_iterator tracks_end() const;
  podio::RelationRange<extension::Track> getTracks() const;
  std::size_t subdetectorHitNumbers_size() const;
  std::int32_t getSubdetectorHitNumbers(std::size_t) const;
  std::vector<std::int32_t>::const_iterator subdetectorHitNumbers_begin() const;
  std::vector<std::int32_t>::const_iterator subdetectorHitNumbers_end() const;
  podio::RelationRange<std::int32_t> getSubdetectorHitNumbers() const;
  std::size_t trackStates_size() const;
  edm4hep::TrackState getTrackStates(std::size_t) const;
  std::vector<edm4hep::TrackState>::const_iterator trackStates_begin() const;
  std::vector<edm4hep::TrackState>::const_iterator trackStates_end() const;
  podio::RelationRange<edm4hep::TrackState> getTrackStates() const;
  std::size_t dxQuantities_size() const;
  edm4hep::Quantity getDxQuantities(std::size_t) const;
  std::vector<edm4hep::Quantity>::const_iterator dxQuantities_begin() const;
  std::vector<edm4hep::Quantity>::const_iterator dxQuantities_end() const;
  podio::RelationRange<edm4hep::Quantity> getDxQuantities() const;


  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from TrackObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<TrackObj>{nullptr}; }

  constexpr bool operator==(const Track& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableTrack& other) const;

  constexpr bool operator!=(const Track& other) const { return !(*this == other); }
  bool operator!=(const MutableTrack& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const Track& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<Track>;

  friend void swap(Track& a, Track& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing TrackObj
  explicit Track(podio::utils::MaybeSharedPtr<TrackObj> obj);
  Track(TrackObj* obj);

  podio::utils::MaybeSharedPtr<TrackObj> m_obj{new TrackObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const Track& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const Track& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::Track> {
  std::size_t operator()(const extension::Track& obj) const {
    return std::hash<extension::TrackObj*>{}(obj.m_obj.get());
  }
};


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::Track::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::Track::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_Track_CollectionData_H
#define EXTENSION_Track_CollectionData_H

// datamodel specific includes
#include "extension/TrackData.h"
#include "extension/TrackObj.h"
#include "extension/Track.h"
#include "extension/TrackerHit.h"

// schema evolution specific includes

// podio specific includes
#include "podio/CollectionBuffers.h"
#include "podio/ICollectionProvider.h"

#include <deque>
#include <memory>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


using TrackObjPointerContainer = std::deque<TrackObj*>;
using TrackDataContainer = std::vector<TrackData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class TrackCollectionData {
public:
  /**
   * The Objs of this collection
   */
  TrackObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  TrackCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  TrackCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  TrackCollectionData(const TrackCollectionData&) = delete;
  TrackCollectionData& operator=(const TrackCollectionData&) = delete;
  TrackCollectionData(TrackCollectionData&& other) = default;
  TrackCollectionData& operator=(TrackCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~TrackCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();

  void createRelations(TrackObj* obj);

  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations
  podio::UVecPtr<extension::TrackerHit> m_rel_trackerHits{std::make_unique<std::vector<extension::TrackerHit>>()}; ///< Relation buffer for read / write
  std::vector<podio::UVecPtr<extension::TrackerHit>> m_rel_trackerHits_tmp{}; ///< Relation buffer for internal book-keeping
  podio::UVecPtr<extension::Track> m_rel_tracks{std::make_unique<std::vector<extension::Track>>()}; ///< Relation buffer for read / write
  std::vector<podio::UVecPtr<extension::Track>> m_rel_tracks_tmp{}; ///< Relation buffer for internal book-keeping

  // members to handle vector members
  podio::UVecPtr<std::int32_t> m_vec_subdetectorHitNumbers{nullptr}; /// combined vector of all objects in collection
  std::vector<podio::UVecPtr<std::int32_t>> m_vecs_subdetectorHitNumbers{}; /// pointers to individual member vectors
  podio::UVecPtr<edm4hep::TrackState> m_vec_trackStates{nullptr}; /// combined vector of all objects in collection
  std::vector<podio::UVecPtr<edm4hep::TrackState>> m_vecs_trackStates{}; /// pointers to individual member vectors
  podio::UVecPtr<edm4hep::Quantity> m_vec_dxQuantities{nullptr}; /// combined vector of all objects in collection
  std::vector<podio::UVecPtr<edm4hep::Quantity>> m_vecs_dxQuantities{}; /// pointers to individual member vectors

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<TrackDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_TrackDATA_H
#define EXTENSION_TrackDATA_H

#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class TrackData
 *  Reconstructed track
 *  @author: EDM4hep authors
 */
class TrackData {
public:
  std::int32_t type{}; ///< flagword that defines the type of track.Bits 16-31 are used internally
  float chi2{}; ///< Chi^2 of the track fit
  std::int32_t ndf{}; ///< number of degrees of freedom of the track fit
  float dEdx{}; ///< dEdx of the track
  float dEdxError{}; ///< error of dEdx
  float radiusOfInnermostHit{}; ///< radius of the innermost hit that has been used in the track fit

  unsigned int subdetectorHitNumbers_begin{};
  unsigned int subdetectorHitNumbers_end{};
  unsigned int trackStates_begin{};
  unsigned int trackStates_end{};
  unsigned int dxQuantities_begin{};
  unsigned int dxQuantities_end{};
  unsigned int trackerHits_begin{};
  unsigned int trackerHits_end{};
  unsigned int tracks_begin{};
  unsigned int tracks_end{};
};


namespace v1 {
using TrackData = extension::TrackData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_TrackOBJ_H
#define EXTENSION_TrackOBJ_H

// data model specific includes
#include "extension/TrackData.h"
#include "extension/TrackerHit.h"
#include "edm4hep/Quantity.h"
#include "edm4hep/TrackState.h"
#include <vector>

#include "podio/ObjectID.h"
#include <vector>



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class Track;

class TrackObj {
public:
  /// constructor
  TrackObj();
  /// copy constructor (does a deep-copy of relation containers)
  TrackObj(const TrackObj&);
  /// constructor from ObjectID and TrackData
  /// does not initialize the internal relation containers
  TrackObj(const podio::ObjectID& id, const TrackData& data);
  /// No assignment operator
  TrackObj& operator=(const TrackObj&) = delete;
  virtual ~TrackObj();

public:
  podio::ObjectID id{};
  TrackData data;
  std::vector<extension::TrackerHit>* m_trackerHits{nullptr};
  std::vector<extension::Track>* m_tracks{nullptr};
  std::vector<std::int32_t>* m_subdetectorHitNumbers{nullptr};
  std::vector<edm4hep::TrackState>* m_trackStates{nullptr};
  std::vector<edm4hep::Quantity>* m_dxQuantities{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_TrackerHit_H
#define EXTENSION_TrackerHit_H

#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/DriftChamberDigiCollection.h"

#include "podio/ObjectID.h"
#include "podio/utilities/TypeHelpers.h"
#include "podio/detail/OrderKey.h"

#include <memory>
#include <ostream>
#include <stdexcept>

namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class TrackerHit
 *  Tracker hit interface class
 *  @author: Thomas Madlener, DESY
 */
class TrackerHit {
public:
  /// type alias containing all the types this interface should work for in a
  /// tuple
  using interfaced_types = std::tuple<::edm4hep::TrackerHitPlane, ::edm4hep::TrackerHit3D, ::extension::SenseWireHit, ::extension::DriftChamberDigi>;

  /// Mark interfaces as always acting as default datatypes.
  ///
  /// @note We only detect the member typedef, not whether it is actually usable
  /// in any meaningful way
  using mutable_type = podio::det::nonesuch;

private:
  /// type alias containing all the mutable types that can be used to initialize
  /// this interface
  using InterfacedMutableTypes = podio::detail::TupleOfMutableTypes<interfaced_types>;

  /// template variable for determining whether type T is a valid interface type
  template<typename T>
  constexpr static bool isInterfacedType = podio::detail::isInTuple<T, interfaced_types>;

public:
  /// template variable for determining whether type T can be used to initialize
  /// this interface
  template<typename T>
  constexpr static bool isInitializableFrom = isInterfacedType<T> || podio::detail::isInTuple<T, InterfacedMutableTypes>;

private:
  struct Concept {
    virtual ~Concept() = default;
    virtual std::unique_ptr<Concept> clone() const = 0;
    virtual void print(std::ostream&) const = 0;

    virtual podio::ObjectID getObjectID() const = 0;
    virtual bool isAvailable() const = 0;
    virtual void unlink() = 0;
    virtual std::uint64_t getCellID() const = 0;

    virtual std::int32_t getType() const = 0;

    virtual std::int32_t getQuality() const = 0;

    virtual float getTime() const = 0;

    virtual float getEDep() const = 0;

    virtual float getEDepError() const = 0;

    virtual const edm4hep::Vector3d& getPosition() const = 0;


    virtual const std::type_info& typeInfo() const = 0;
    virtual bool equal(const Concept* rhs) const = 0;
    virtual podio::detail::OrderKey objOrderKey() const = 0;
    virtual size_t objHash() const = 0;
  };

  template<typename ValueT>
  struct Model final : Concept {
    ~Model() override = default;
    Model(const ValueT& value) : m_value(value) {}

    std::unique_ptr<Concept> clone() const final {
      return std::make_unique<Model<ValueT>>(m_value);
    }

    void print(std::ostream& os) const final {
      os << m_value;
    }

    void unlink() final { m_value.unlink(); }
    bool isAvailable() const final { return m_value.isAvailable(); }
    podio::ObjectID getObjectID() const final { return m_value.getObjectID(); }

    const std::type_info& typeInfo() const final { return typeid(ValueT); }

    bool equal(const Concept* rhs) const final {
      if (typeInfo() == rhs->typeInfo()) {
        return m_value == static_cast<const Model<ValueT>*>(rhs)->m_value;
      }
      return false;
    }

    podio::detail::OrderKey objOrderKey() const final {
      return podio::detail::getOrderKey(m_value);
    }

    size_t objHash() const final {return std::hash<ValueT>{}(m_value); }

    std::uint64_t getCellID() const final { return m_value.getCellID(); }

    std::int32_t getType() const final { return m_value.getType(); }

    std::int32_t getQuality() const final { return m_value.getQuality(); }

    float getTime() const final { return m_value.getTime(); }

    float getEDep() const final { return m_value.getEDep(); }

    float getEDepError() const final { return m_value.getEDepError(); }

    const edm4hep::Vector3d& getPosition() const final { return m_value.getPosition(); }



    ValueT m_value{};
  };

 std::unique_ptr<Concept> m_self{nullptr};

public:
  // TrackerHit can only be initialized with one of the following types (and their Mutable counter parts): ::edm4hep::TrackerHitPlane, ::edm4hep::TrackerHit3D, ::extension::SenseWireHit, ::extension::DriftChamberDigi
  template<typename ValueT>
  requires isInitializableFrom<ValueT>
  TrackerHit(const ValueT& value) :
    m_self(std::make_unique<Model<podio::detail::GetDefaultHandleType<ValueT>>>(value)) {
  }

  TrackerHit(const TrackerHit& other) :
    m_self(other.m_self->clone()) {}
  TrackerHit& operator=(const TrackerHit& other) {
    TrackerHit tmp{other};
    std::swap(tmp.m_self, this->m_self);
    return *this;
  }

  ~TrackerHit() = default;
  TrackerHit(TrackerHit&&) = default;
  TrackerHit& operator=(TrackerHit&&) = default;

  /// Create an empty handle
  static TrackerHit makeEmpty() {
    // We simply chose the first type of the interfaced types here to initialize
    // an empty handle
    return ::edm4hep::TrackerHitPlane::makeEmpty();
  }

  static constexpr std::string_view typeName = "extension::TrackerHit";

  /// check whether the object is actually available
  bool isAvailable() const { return m_self->isAvailable(); }
  /// disconnect from the underlying value
  void unlink() { m_self->unlink(); }

  podio::ObjectID id() const { return getObjectID(); }
  podio::ObjectID getObjectID() const { return m_self->getObjectID(); }

  /// Check if the object currently holds a value of the requested type
  template<typename T>
  bool isA() const {
    static_assert(isInterfacedType<T>, "TrackerHit can only ever be one of the following types: ::edm4hep::TrackerHitPlane, ::edm4hep::TrackerHit3D, ::extension::SenseWireHit, ::extension::DriftChamberDigi");
    return typeid(T) == m_self->typeInfo();
  }

  /// Get the contained value as the concrete type it was put in. This will
  /// throw a std::runtime_error if T is not the type of the currently held
  /// value. Use isA to check beforehand if necessary
  template<typename T>
  T as() const {
    if (!isA<T>()) {
      throw std::runtime_error("Cannot get value as object currently holds another type");
    }
    // We can safely cast here since we check types before
    return static_cast<Model<T>*>(m_self.get())->m_value;
  }

  friend bool operator==(const TrackerHit& lhs, const TrackerHit& rhs) {
    return lhs.m_self->equal(rhs.m_self.get());
  }

  friend bool operator!=(const TrackerHit& lhs, const TrackerHit& rhs) {
    return !(lhs == rhs);
  }

  friend bool operator<(const TrackerHit& lhs, const TrackerHit& rhs) {
    return lhs.m_self->objOrderKey() < rhs.m_self->objOrderKey();
  }

  /// Access the ID of the sensor that created this hit
  std::uint64_t getCellID() const { return m_self->getCellID(); }

  /// Access the type of the raw data hit
  std::int32_t getType() const { return m_self->getType(); }

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const { return m_self->getQuality(); }

  /// Access the time of the hit [ns]
  float getTime() const { return m_self->getTime(); }

  /// Access the energy deposited on the hit [GeV]
  float getEDep() const { return m_self->getEDep(); }

  /// Access the error measured on eDep [GeV]
  float getEDepError() const { return m_self->getEDepError(); }

  /// Access the hit position [mm]
  const edm4hep::Vector3d& getPosition() const { return m_self->getPosition(); }



  friend std::ostream& operator<<(std::ostream& os, const TrackerHit& value) {
    value.m_self->print(os);
    return os;
  }

  friend std::hash<TrackerHit>;
};

/// @}
} // namespace extension


template<>
struct std::hash<extension::TrackerHit> {
  std::size_t operator()(const extension::TrackerHit& obj) const {
    return obj.m_self->objHash();
  }
};

#endif
// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_TrackerHit_dev_H
#define EXTENSION_TrackerHit_dev_H

#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "extension/SenseWireHitCollection.h"

#include "podio/ObjectID.h"
#include "podio/utilities/TypeHelpers.h"
#include "podio/detail/OrderKey.h"

#include <memory>
#include <ostream>
#include <stdexcept>

namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class TrackerHit_dev
 *  Tracker hit interface class to test the mixing of extension and edm4hep data types
 *  @author: Thomas Madlener, DESY
 */
class TrackerHit_dev {
public:
  /// type alias containing all the types this interface should work for in a
  /// tuple
  using interfaced_types = std::tuple<::edm4hep::TrackerHitPlane, ::edm4hep::TrackerHit3D, ::extension::SenseWireHit>;

  /// Mark interfaces as always acting as default datatypes.
  ///
  /// @note We only detect the member typedef, not whether it is actually usable
  /// in any meaningful way
  using mutable_type = podio::det::nonesuch;

private:
  /// type alias containing all the mutable types that can be used to initialize
  /// this interface
  using InterfacedMutableTypes = podio::detail::TupleOfMutableTypes<interfaced_types>;

  /// template variable for determining whether type T is a valid interface type
  template<typename T>
  constexpr static bool isInterfacedType = podio::detail::isInTuple<T, interfaced_types>;

public:
  /// template variable for determining whether type T can be used to initialize
  /// this interface
  template<typename T>
  constexpr static bool isInitializableFrom = isInterfacedType<T> || podio::detail::isInTuple<T, InterfacedMutableTypes>;

private:
  struct Concept {
    virtual ~Concept() = default;
    virtual std::unique_ptr<Concept> clone() const = 0;
    virtual void print(std::ostream&) const = 0;

    virtual podio::ObjectID getObjectID() const = 0;
    virtual bool isAvailable() const = 0;
    virtual void unlink() = 0;
    virtual std::uint64_t getCellID() const = 0;

    virtual std::int32_t getType() const = 0;

    virtual std::int32_t getQuality() const = 0;

    virtual float getTime() const = 0;

    virtual float getEDep() const = 0;

    virtual float getEDepError() const = 0;

    virtual const edm4hep::Vector3d& getPosition() const = 0;


    virtual const std::type_info& typeInfo() const = 0;
    virtual bool equal(const Concept* rhs) const = 0;
    virtual podio::detail::OrderKey objOrderKey() const = 0;
    virtual size_t objHash() const = 0;
  };

  template<typename ValueT>
  struct Model final : Concept {
    ~Model() override = default;
    Model(const ValueT& value) : m_value(value) {}

    std::unique_ptr<Concept> clone() const final {
      return std::make_unique<Model<ValueT>>(m_value);
    }

    void print(std::ostream& os) const final {
      os << m_value;
    }

    void unlink() final { m_value.unlink(); }
    bool isAvailable() const final { return m_value.isAvailable(); }
    podio::ObjectID getObjectID() const final { return m_value.getObjectID(); }

    const std::type_info& typeInfo() const final { return typeid(ValueT); }

    bool equal(const Concept* rhs) const final {
      if (typeInfo() == rhs->typeInfo()) {
        return m_value == static_cast<const Model<ValueT>*>(rhs)->m_value;
      }
      return false;
    }

    podio::detail::OrderKey objOrderKey() const final {
      return podio::detail::getOrderKey(m_value);
    }

    size_t objHash() const final {return std::hash<ValueT>{}(m_value); }

    std::uint64_t getCellID() const final { return m_value.getCellID(); }

    std::int32_t getType() const final { return m_value.getType(); }

    std::int32_t getQuality() const final { return m_value.getQuality(); }

    float getTime() const final { return m_value.getTime(); }

    float getEDep() const final { return m_value.getEDep(); }

    float getEDepError() const final { return m_value.getEDepError(); }

    const edm4hep::Vector3d& getPosition() const final { return m_value.getPosition(); }



    ValueT m_value{};
  };

 std::unique_ptr<Concept> m_self{nullptr};

public:
  // TrackerHit_dev can only be initialized with one of the following types (and their Mutable counter parts): ::edm4hep::TrackerHitPlane, ::edm4hep::TrackerHit3D, ::extension::SenseWireHit
  template<typename ValueT>
  requires isInitializableFrom<ValueT>
  TrackerHit_dev(const ValueT& value) :
    m_self(std::make_unique<Model<podio::detail::GetDefaultHandleType<ValueT>>>(value)) {
  }

  TrackerHit_dev(const TrackerHit_dev& other) :
    m_self(other.m_self->clone()) {}
  TrackerHit_dev& operator=(const TrackerHit_dev& other) {
    TrackerHit_dev tmp{other};
    std::swap(tmp.m_self, this->m_self);
    return *this;
  }

  ~TrackerHit_dev() = default;
  TrackerHit_dev(TrackerHit_dev&&) = default;
  TrackerHit_dev& operator=(TrackerHit_dev&&) = default;

  /// Create an empty handle
  static TrackerHit_dev makeEmpty() {
    // We simply chose the first type of the interfaced types here to initialize
    // an empty handle
    return ::edm4hep::TrackerHitPlane::makeEmpty();
  }

  static constexpr std::string_view typeName = "extension::TrackerHit_dev";

  /// check whether the object is actually available
  bool isAvailable() const { return m_self->isAvailable(); }
  /// disconnect from the underlying value
  void unlink() { m_self->unlink(); }

  podio::ObjectID id() const { return getObjectID(); }
  podio::ObjectID getObjectID() const { return m_self->getObjectID(); }

  /// Check if the object currently holds a value of the requested type
  template<typename T>
  bool isA() const {
    static_assert(isInterfacedType<T>, "TrackerHit_dev can only ever be one of the following types: ::edm4hep::TrackerHitPlane, ::edm4hep::TrackerHit3D, ::extension::SenseWireHit");
    return typeid(T) == m_self->typeInfo();
  }

  /// Get the contained value as the concrete type it was put in. This will
  /// throw a std::runtime_error if T is not the type of the currently held
  /// value. Use isA to check beforehand if necessary
  template<typename T>
  T as() const {
    if (!isA<T>()) {
      throw std::runtime_error("Cannot get value as object currently holds another type");
    }
    // We can safely cast here since we check types before
    return static_cast<Model<T>*>(m_self.get())->m_value;
  }

  friend bool operator==(const TrackerHit_dev& lhs, const TrackerHit_dev& rhs) {
    return lhs.m_self->equal(rhs.m_self.get());
  }

  friend bool operator!=(const TrackerHit_dev& lhs, const TrackerHit_dev& rhs) {
    return !(lhs == rhs);
  }

  friend bool operator<(const TrackerHit_dev& lhs, const TrackerHit_dev& rhs) {
    return lhs.m_self->objOrderKey() < rhs.m_self->objOrderKey();
  }

  /// Access the ID of the sensor that created this hit
  std::uint64_t getCellID() const { return m_self->getCellID(); }

  /// Access the type of the raw data hit
  std::int32_t getType() const { return m_self->getType(); }

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const { return m_self->getQuality(); }

  /// Access the time of the hit [ns]
  float getTime() const { return m_self->getTime(); }

  /// Access the energy deposited on the hit [GeV]
  float getEDep() const { return m_self->getEDep(); }

  /// Access the error measured on eDep [GeV]
  float getEDepError() const { return m_self->getEDepError(); }

  /// Access the hit position [mm]
  const edm4hep::Vector3d& getPosition() const { return m_self->getPosition(); }



  friend std::ostream& operator<<(std::ostream& os, const TrackerHit_dev& value) {
    value.m_self->print(os);
    return os;
  }

  friend std::hash<TrackerHit_dev>;
};

/// @}
} // namespace extension


template<>
struct std::hash<extension::TrackerHit_dev> {
  std::size_t operator()(const extension::TrackerHit_dev& obj) const {
    return obj.m_self->objHash();
  }
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"extension::DriftChamberDigi", payloadCode, "@",
"extension::DriftChamberDigiCollection", payloadCode, "@",
"extension::DriftChamberDigiData", payloadCode, "@",
"extension::DriftChamberDigiLocal", payloadCode, "@",
"extension::DriftChamberDigiLocalCollection", payloadCode, "@",
"extension::DriftChamberDigiLocalData", payloadCode, "@",
"extension::DriftChamberDigiV2", payloadCode, "@",
"extension::DriftChamberDigiV2Collection", payloadCode, "@",
"extension::DriftChamberDigiV2Data", payloadCode, "@",
"extension::MCRecoDriftChamberDigiAssociation", payloadCode, "@",
"extension::MCRecoDriftChamberDigiAssociationCollection", payloadCode, "@",
"extension::MCRecoDriftChamberDigiAssociationData", payloadCode, "@",
"extension::MCRecoDriftChamberDigiV2Association", payloadCode, "@",
"extension::MCRecoDriftChamberDigiV2AssociationCollection", payloadCode, "@",
"extension::MCRecoDriftChamberDigiV2AssociationData", payloadCode, "@",
"extension::MutableDriftChamberDigi", payloadCode, "@",
"extension::MutableDriftChamberDigiLocal", payloadCode, "@",
"extension::MutableDriftChamberDigiV2", payloadCode, "@",
"extension::MutableMCRecoDriftChamberDigiAssociation", payloadCode, "@",
"extension::MutableMCRecoDriftChamberDigiV2Association", payloadCode, "@",
"extension::MutableSenseWireHit", payloadCode, "@",
"extension::MutableSenseWireHitSimTrackerHitLink", payloadCode, "@",
"extension::MutableTrack", payloadCode, "@",
"extension::SenseWireHit", payloadCode, "@",
"extension::SenseWireHitCollection", payloadCode, "@",
"extension::SenseWireHitData", payloadCode, "@",
"extension::SenseWireHitSimTrackerHitLink", payloadCode, "@",
"extension::SenseWireHitSimTrackerHitLinkCollection", payloadCode, "@",
"extension::SenseWireHitSimTrackerHitLinkData", payloadCode, "@",
"extension::Track", payloadCode, "@",
"extension::TrackCollection", payloadCode, "@",
"extension::TrackData", payloadCode, "@",
"extension::v1::DriftChamberDigiData", payloadCode, "@",
"extension::v1::DriftChamberDigiLocalData", payloadCode, "@",
"extension::v1::DriftChamberDigiV2Data", payloadCode, "@",
"extension::v1::MCRecoDriftChamberDigiAssociationData", payloadCode, "@",
"extension::v1::MCRecoDriftChamberDigiV2AssociationData", payloadCode, "@",
"extension::v1::SenseWireHitData", payloadCode, "@",
"extension::v1::SenseWireHitSimTrackerHitLinkData", payloadCode, "@",
"extension::v1::TrackData", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libextensionDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libextensionDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libextensionDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libextensionDict() {
  TriggerDictionaryInitialization_libextensionDict_Impl();
}
