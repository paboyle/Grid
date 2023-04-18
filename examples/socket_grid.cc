#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>
#include <stdio.h>
#include <err.h>
#include <fcntl.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

static int sock;
static const char *sock_path_fmt = "/tmp/GridUnixSocket.%d";
static char sock_path[256];

class UnixSockets {
public:
  static void Open(int rank)
  {
    int errnum;

    sock = socket(AF_UNIX, SOCK_DGRAM, 0);  assert(sock>0);
    printf("allocated socket %d\n",sock);

    struct sockaddr_un sa_un = { 0 };
    sa_un.sun_family = AF_UNIX;
    snprintf(sa_un.sun_path, sizeof(sa_un.sun_path),sock_path_fmt,rank);
    unlink(sa_un.sun_path);
    if (bind(sock, (struct sockaddr *)&sa_un, sizeof(sa_un))) {
      perror("bind failure");
      exit(EXIT_FAILURE);
    }
    printf("bound socket %d to %s\n",sock,sa_un.sun_path);
  }

  static int RecvFileDescriptor(void)
  {
    int n;
    int fd;
    char buf[1];
    struct iovec iov;
    struct msghdr msg;
    struct cmsghdr *cmsg;
    char cms[CMSG_SPACE(sizeof(int))];

    iov.iov_base = buf;
    iov.iov_len = 1;

    memset(&msg, 0, sizeof msg);
    msg.msg_name = 0;
    msg.msg_namelen = 0;
    msg.msg_iov = &iov;
    msg.msg_iovlen = 1;

    msg.msg_control = (caddr_t)cms;
    msg.msg_controllen = sizeof cms;

    if((n=recvmsg(sock, &msg, 0)) < 0) {
      perror("recvmsg failed");
      return -1;
    }
    if(n == 0){
      perror("recvmsg returned 0");
      return -1;
    }
    cmsg = CMSG_FIRSTHDR(&msg);
    memmove(&fd, CMSG_DATA(cmsg), sizeof(int));
    printf("received fd %d from socket %d\n",fd,sock);
    return fd;
  }

  static void SendFileDescriptor(int fildes,int xmit_to_rank)
  {
    struct msghdr msg;
    struct iovec iov;
    struct cmsghdr *cmsg = NULL;
    char ctrl[CMSG_SPACE(sizeof(int))];
    char data = ' ';

    memset(&msg, 0, sizeof(struct msghdr));
    memset(ctrl, 0, CMSG_SPACE(sizeof(int)));
    iov.iov_base = &data;
    iov.iov_len = sizeof(data);
    
    sprintf(sock_path,sock_path_fmt,xmit_to_rank);
    printf("sending FD %d over socket %d to rank %d AF_UNIX path %s\n",fildes,sock,xmit_to_rank,sock_path);fflush(stdout);
    
    struct sockaddr_un sa_un = { 0 };
    sa_un.sun_family = AF_UNIX;
    snprintf(sa_un.sun_path, sizeof(sa_un.sun_path),sock_path_fmt,xmit_to_rank);

    msg.msg_name = (void *)&sa_un;
    msg.msg_namelen = sizeof(sa_un);
    msg.msg_iov = &iov;
    msg.msg_iovlen = 1;
    msg.msg_controllen =  CMSG_SPACE(sizeof(int));
    msg.msg_control = ctrl;

    cmsg = CMSG_FIRSTHDR(&msg);
    cmsg->cmsg_level = SOL_SOCKET;
    cmsg->cmsg_type = SCM_RIGHTS;
    cmsg->cmsg_len = CMSG_LEN(sizeof(int));

    *((int *) CMSG_DATA(cmsg)) = fildes;

    if ( sendmsg(sock, &msg, 0) == -1 ) perror("sendmsg failed");
  };
};

int main(int argc, char **argv)
{
  int me = fork()?0:1;
  
  UnixSockets::Open(me);
  
  // need MPI barrier
  sleep(10);
  const char * message = "Hello, World\n";
  if( me ) {
    int fd = open("foo",O_RDWR|O_CREAT,0666);
    if ( fd < 0 ) {
      perror("failed to open file");
      exit(EXIT_FAILURE);
    }
    // rank 1 sends ot rank 0
    UnixSockets::SendFileDescriptor(fd,0);
    close(fd);
  } else {
    // rank 0 sends receives frmo rank 1
    int fd = UnixSockets::RecvFileDescriptor();
    write(fd,(const void *)message,strlen(message));
    close(fd);
  }
}
