#pragma once
#define SECURITY_WIN32
#include <security.h>
#include <Windows.h>

KSECDDDECLSPEC
SECURITY_STATUS SEC_ENTRY
secur32$AcquireCredentialsHandleA(
    _In_opt_  LPSTR pszPrincipal,                 // Name of principal
    _In_      LPSTR pszPackage,                   // Name of package
    _In_      unsigned long fCredentialUse,       // Flags indicating use
    _In_opt_  void * pvLogonId,                   // Pointer to logon ID
    _In_opt_  void * pAuthData,                   // Package specific data
    _In_opt_  SEC_GET_KEY_FN pGetKeyFn,           // Pointer to GetKey() func
    _In_opt_  void * pvGetKeyArgument,            // Value to pass to GetKey()
    _Out_     PCredHandle phCredential,           // (out) Cred Handle
    _Out_opt_ PTimeStamp ptsExpiry                // (out) Lifetime (optional)
    );

KSECDDDECLSPEC
SECURITY_STATUS SEC_ENTRY
secur32$InitializeSecurityContextA(
    _In_opt_    PCredHandle phCredential,               // Cred to base context
    _In_opt_    PCtxtHandle phContext,                  // Existing context (OPT)
    _In_opt_    SEC_CHAR * pszTargetName,       // Name of target
    _In_        unsigned long fContextReq,              // Context Requirements
    _In_        unsigned long Reserved1,                // Reserved, MBZ
    _In_        unsigned long TargetDataRep,            // Data rep of target
    _In_opt_    PSecBufferDesc pInput,                  // Input Buffers
    _In_        unsigned long Reserved2,                // Reserved, MBZ
    _Inout_opt_ PCtxtHandle phNewContext,               // (out) New Context handle
    _Inout_opt_ PSecBufferDesc pOutput,                 // (inout) Output Buffers
    _Out_       unsigned long * pfContextAttr,  // (out) Context attrs
    _Out_opt_   PTimeStamp ptsExpiry                    // (out) Life span (OPT)
    );

KSECDDDECLSPEC
SECURITY_STATUS SEC_ENTRY
secur32$AcceptSecurityContext(
    _In_opt_    PCredHandle phCredential,               // Cred to base context
    _In_opt_    PCtxtHandle phContext,                  // Existing context (OPT)
    _In_opt_    PSecBufferDesc pInput,                  // Input buffer
    _In_        unsigned long fContextReq,              // Context Requirements
    _In_        unsigned long TargetDataRep,            // Target Data Rep
    _Inout_opt_ PCtxtHandle phNewContext,               // (out) New context handle
    _Inout_opt_ PSecBufferDesc pOutput,                 // (inout) Output buffers
    _Out_       unsigned long * pfContextAttr,  // (out) Context attributes
    _Out_opt_   PTimeStamp ptsExpiry                    // (out) Life span (OPT)
    );

KSECDDDECLSPEC
SECURITY_STATUS SEC_ENTRY
secur32$QuerySecurityContextToken(
	_In_  PCtxtHandle phContext,
	_Out_ void **Token
);

KSECDDDECLSPEC
SECURITY_STATUS SEC_ENTRY
secur32$FreeCredentialsHandle(
    _In_ PCredHandle phCredential            // Handle to free
    );

KSECDDDECLSPEC
SECURITY_STATUS SEC_ENTRY
secur32$FreeContextBuffer(
    _Inout_ PVOID pvContextBuffer      // buffer to free
    );

KSECDDDECLSPEC
SECURITY_STATUS SEC_ENTRY
secur32$DeleteSecurityContext(
    _In_ PCtxtHandle phContext               // Context to delete
    );
