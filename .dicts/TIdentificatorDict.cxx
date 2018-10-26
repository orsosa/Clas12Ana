// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdictsdITIdentificatorDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TIdentificator.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *TIdentificator_Dictionary();
   static void TIdentificator_TClassManip(TClass*);
   static void *new_TIdentificator(void *p = 0);
   static void *newArray_TIdentificator(Long_t size, void *p);
   static void delete_TIdentificator(void *p);
   static void deleteArray_TIdentificator(void *p);
   static void destruct_TIdentificator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TIdentificator*)
   {
      ::TIdentificator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::TIdentificator));
      static ::ROOT::TGenericClassInfo 
         instance("TIdentificator", "TIdentificator.h", 34,
                  typeid(::TIdentificator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TIdentificator_Dictionary, isa_proxy, 0,
                  sizeof(::TIdentificator) );
      instance.SetNew(&new_TIdentificator);
      instance.SetNewArray(&newArray_TIdentificator);
      instance.SetDelete(&delete_TIdentificator);
      instance.SetDeleteArray(&deleteArray_TIdentificator);
      instance.SetDestructor(&destruct_TIdentificator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TIdentificator*)
   {
      return GenerateInitInstanceLocal((::TIdentificator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TIdentificator*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TIdentificator_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TIdentificator*)0x0)->GetClass();
      TIdentificator_TClassManip(theClass);
   return theClass;
   }

   static void TIdentificator_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_TIdentificator(void *p) {
      return  p ? new(p) ::TIdentificator : new ::TIdentificator;
   }
   static void *newArray_TIdentificator(Long_t nElements, void *p) {
      return p ? new(p) ::TIdentificator[nElements] : new ::TIdentificator[nElements];
   }
   // Wrapper around operator delete
   static void delete_TIdentificator(void *p) {
      delete ((::TIdentificator*)p);
   }
   static void deleteArray_TIdentificator(void *p) {
      delete [] ((::TIdentificator*)p);
   }
   static void destruct_TIdentificator(void *p) {
      typedef ::TIdentificator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TIdentificator

namespace {
  void TriggerDictionaryInitialization_TIdentificatorDict_Impl() {
    static const char* headers[] = {
"TIdentificator.h",
0
    };
    static const char* includePaths[] = {
"/home/orsosa/ClasTool/include",
"/opt/root6/include",
"/home/orsosa/Analyser/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TIdentificatorDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TIdentificator.h")))  TIdentificator;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TIdentificatorDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TIdentificator.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TIdentificator", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TIdentificatorDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TIdentificatorDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TIdentificatorDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TIdentificatorDict() {
  TriggerDictionaryInitialization_TIdentificatorDict_Impl();
}
