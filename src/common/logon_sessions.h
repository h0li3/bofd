#pragma once
#include "beacon.h"
#include <ntsecapi.h>

DECLSPEC_IMPORT NTSTATUS __stdcall secur32$LsaEnumerateLogonSessions(
	_Out_ PULONG LogonSessionCount,
	_Out_ PLUID  *LogonSessionList
);

DECLSPEC_IMPORT NTSTATUS __stdcall secur32$LsaGetLogonSessionData(
	_Out_ PLUID                        LogonId,
	_Out_ PSECURITY_LOGON_SESSION_DATA *ppLogonSessionData
);

DECLSPEC_IMPORT NTSTATUS __stdcall secur32$LsaFreeReturnBuffer(
	_In_ PVOID Buffer
);

DECLSPEC_IMPORT ULONG __stdcall advapi32$LsaNtStatusToWinError(
	_In_ NTSTATUS Status
);

typedef BOOL (*ListLogonSessionsCallback)(SECURITY_LOGON_SESSION_DATA *, void* user_data);

void list_logon_sessions(ListLogonSessionsCallback callback, void* user_data)
{
	NTSTATUS status;
	ULONG session_count = 0;
	LUID *session_list = NULL;
	SECURITY_LOGON_SESSION_DATA* session_data;

	status = secur32$LsaEnumerateLogonSessions(&session_count, &session_list);
	if (status != 0) {
		BeaconPrintf(CALLBACK_ERROR, "could not enumerate sessions: 0x%08lx", advapi32$LsaNtStatusToWinError(status));
		return;
	}

	for (ULONG i = 0; i < session_count; ++i) {
		status = secur32$LsaGetLogonSessionData(&session_list[i], &session_data);
		if (status != 0) {
			//BeaconPrintf(CALLBACK_ERROR, "could not get session data: 0x%08lx", advapi32$LsaNtStatusToWinError(status));
			continue;
		}
		BOOL goon = callback(session_data, user_data);
		secur32$LsaFreeReturnBuffer(session_data);
		if (!goon) break;
		session_data = 0;
	}

	secur32$LsaFreeReturnBuffer(session_list);
}
