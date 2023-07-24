#include <WinSock2.h>
#include "common/beacon.h"
#include <string.h>
#include <stdlib.h>

int parse_address(sockaddr_in& addr, const char* host, int port)
{
	auto ent = gethostbyname(host);
	if (ent->h_length == 0) {
		return -1;
	}
	addr.sin_family = AF_INET;
	addr.sin_addr.S_un.S_addr = *(ULONG*)ent->h_addr_list[0];
	addr.sin_port = htons(port);
	return 0;
}

extern "C" void go(char* arg, int alen)
{
	int ret;
	char* host;
	int port;
	sockaddr_in addr;
	char* p = nullptr;

	for (int i = 0; i < alen; ++i) {
		if (arg[i] == ':') {
			p = arg + i;
			break;
		}
	}

	if (p == nullptr || p + 1 == arg + alen) {
        bof::errorf("illegal target: %s", arg);
		return;
	}
	*p = 0;
	host = arg;
	port = atoi(p+1);

	WSADATA wd;
	ret = WSAStartup(MAKEWORD(1, 1), &wd);
	if (ret == SOCKET_ERROR) {
        bof::errorf("failed to startup socket: %d", ret);
		return;
	}

	SOCKET sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
	if (sock == SOCKET_ERROR) {
        bof::errorf("failed to create socket: %d", sock);
		return;
	}

	ret = parse_address(addr, host, port);
	if (ret != 0) {
        bof::errorf("failed to resolve host: %s", host);
		closesocket(sock);
		return;
	}

	unsigned long ul = 1;
	ret = ioctlsocket(sock, FIONBIO, (unsigned long*)&ul);
	if (ret == SOCKET_ERROR) return;

	fd_set r;
	FD_ZERO(&r);
	FD_SET(sock, &r);

	ret = connect(sock, (const sockaddr*)&addr, sizeof(addr));
	timeval timeout {1, 0};
	ret = select(0, 0, &r, 0, &timeout);
	if (ret <= 0)
	{
		closesocket(sock);
		return;
	}

    bof::printf("%s:%d    open", host, port);

	closesocket(sock);
}
