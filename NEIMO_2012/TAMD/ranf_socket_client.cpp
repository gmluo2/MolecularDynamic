#include "project.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if _SYS_ == _WINDOWS_SYS_
#include <Windows.h>
#include "ranlib.h"
int connect_ranf_server(const char *host) {return -1;}
float get_ranf() {return (float)ranf();};
void disconnect_ranf_server() {return;}

#elif _SYS_ == _LINUX_SYS_
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 

void error(const char *msg)
{
  perror(msg);
  exit(0);
}

int sockfd = -1;
extern void show_msg(char *buff, bool endline = true);

void local_show_msg(char *buff, bool endline = true) {
  //printf("%s", buff); if (endline) printf("\n");
  return show_msg(buff, endline);
}

int connect_ranf_server(const char *host) {
  //int sockfd;
  int  portno, n;
  struct sockaddr_in serv_addr;
  struct hostent *server;
  char buffer[256];
  //if (argc < 3) {
  //  fprintf(stderr,"usage %s hostname port\n", argv[0]);
  //  exit(0);
  //}
  //portno = atoi(argv[2]);
  portno=51717; // port specifically for universal ranf service
  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd < 0) error("ERROR opening socket");
  server = gethostbyname(host);
  if (server == NULL) {
    fprintf(stderr,"ERROR, no such host\n");
    //exit(0);
	return -1;
  }
  bzero((char *) &serv_addr, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  bcopy((char *)server->h_addr,
        (char *)&serv_addr.sin_addr.s_addr,
        server->h_length);
  serv_addr.sin_port = htons(portno);

  int itry = 0, ntry = 20, nwait = 2; // wait 2s between each try
  while (itry < ntry) {
    if (connect(sockfd,(struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0) {itry++; sleep(nwait); continue;} //error("ERROR connecting");
    else return 1;
  }
  sockfd = -1; return -1; // failure to setup the connection
}

float get_ranf() {
  if (sockfd < 0) return 0;
  char buff[512] = "\0", obuff[256] = "\0";
  int n;
  float f = 0;
  strcpy(obuff, "REQUEST_RANF");
  if (send(sockfd, obuff, strlen(obuff) + 1, MSG_OOB) < 0) return 0;
  //printf("send out request\n");
  bzero(buff,512);
  n = recv(sockfd, buff, 510, MSG_WAITALL);
  if (n < 0) local_show_msg("failure to get universal random float "); //error("ERROR reading from socket");
  //printf("%s\n",buffer);
  //close(sockfd);
  else if (sscanf(buff, "%f", &f) != 1) {
    local_show_msg("unknown universal random float: ", false); local_show_msg(buff, true);
  }
  //else printf("received %s;  ", buff);
  return f;
}

void disconnect_ranf_server() {
  if (sockfd >= 0) send(sockfd, "REQUEST_CLOSE", 14, MSG_OOB);
  close(sockfd); sockfd = -1;
}

/*
int main(int argc, char *argv[]) {
  int itry = 0, irand = 0, nrands = 20;
  char buff[256] = "\0";
  float f = 0;
  while (itry < 10) {
    if (connect_ranf_server("127.0.0.1") < 0) {show_msg("failure to connect server", true); return 1;}
    irand = 0;
    printf("connection -- %d\n", itry);
    while (irand < nrands) {
      f = get_ranf();
      printf("%d -- %f\n", irand + itry * 100, f);
      //sprintf(buff, "%d -- %f", irand + itry * 100, get_ranf());
      //show_msg(buff, true);
      irand++;
    }
    disconnect_ranf_server();
    itry++;
  }
}
*/

#endif
