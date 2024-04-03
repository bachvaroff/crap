#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <errno.h>
#include <prototypes.h>
#include <sysexits.h>

#define MAX_CHLD	32
#define NCHLD		16

struct chld_status {
	int pid;
	int status;
};

struct chld_status chlds[MAX_CHLD];
int hp, tp, lost_child;
int sync_pipe[2];

void
sigcld(sig)
	int sig;
{
	int pid, status, nb;
	char dummy = 'S';
	
	(void)sig;
	
	pid = wait(&status);
	
	if (pid > 0) {
		if (((hp + 1) % MAX_CHLD) == tp) lost_child++;
		else {
			chlds[hp].pid = pid;
			chlds[hp].status = status;
			hp = (hp + 1) % MAX_CHLD;
		}
		
		if (rdchk(sync_pipe[0]) <= 0) do {
			nb = write(sync_pipe[1], &dummy, sizeof (dummy));
			if (nb < 0) {
				if (errno != EINTR) _exit(EX_OSERR);
			} else if (nb == (sizeof (dummy))) break;
			else _exit(EX_OSERR);
		} while (1);
	} else _exit(EX_OSERR);
	
	if (signal(SIGCLD, sigcld) == SIG_ERR)
		_exit(EX_OSERR);
	
	return;
}

void
prepare_sigcld()
{
	int flags;
	
	hp = 0;
	tp = 0;
	lost_child = 0;
	
	if (pipe(sync_pipe) < 0) {
		perror("pipe");
		exit(EX_OSERR);
	}
	
	if ((flags = fcntl(sync_pipe[0], F_GETFL)) < 0) {
		perror("fcntl");
		exit(EX_OSERR);
	}
	flags |= O_NDELAY;
	if (fcntl(sync_pipe[0], F_SETFL, flags) < 0) {
		perror("fcntl");
		exit(EX_OSERR);
	}

	if ((flags = fcntl(sync_pipe[0], F_GETFL)) < 0) {
		perror("fcntl");
		exit(EX_OSERR);
	}
	flags |= O_NDELAY;
	if (fcntl(sync_pipe[0], F_SETFL, flags) < 0) {
		perror("fcntl");
		exit(EX_OSERR);
	}
	
	if (signal(SIGCLD, sigcld) == SIG_ERR) {
		perror("signal");
		exit(EX_OSERR);
	}
	 	
	return;
}

void
done_sigcld()
{
	(void)signal(SIGCLD, SIG_DFL);
	(void)close(sync_pipe[0]);
	(void)close(sync_pipe[1]);
	return;
}

void
do_something()
{
	int j;
	
	printf("doing something...\n");
	for (j = 0; j < 8; j++)
		(void)sleep(1u);
	
	return;
}

void
do_something_child(task)
	unsigned int task;
{
	printf("child %d task %u\n", getpid(), task);
	(void)sleep(task);
	printf("child %d task %u done\n", getpid(), task);
	
	return;
}

int
main()
{
	int pid, pids[NCHLD], j, term, nb;
	char dummy;
	
	(void)signal(SIGPIPE, SIG_IGN);
	prepare_sigcld();
	
	for (j = 0; j < NCHLD; j++) {
		if ((pid = fork()) < 0) {
			perror("fork");
			pids[j] = pid;
		} else if (pid) {
			printf("child pid (fork) = %d\n", pid);
			pids[j] = pid;
		} else {
			done_sigcld();
			do_something_child((unsigned int)j + ((j > 7) ? 10u : 5u));
			printf("child %d exiting...\n", getpid());
			exit(EX_OK);
		}
	}
	
	for (term = 0; !term; ) {
		do_something();
		
		do {
			if (rdchk(sync_pipe[0]) <= 0) {
				printf("no child to check for...\n");
				break;
			} else do {
				nb = read(sync_pipe[0], &dummy, sizeof (dummy));
				if (nb < 0) {
					if (errno != EINTR) exit(EX_OSERR);
				} else if (nb == (sizeof (dummy))) break;
				else printf("WTF!\n");
			} while (1);
			
			printf("got '%c', tp == %d, hp == %d\n", dummy, tp, hp);
			
			for (; tp != hp; tp = (tp + 1) % MAX_CHLD) {
				printf("\ttp == %d; chlds[tp].pid == %d; chlds[tp].status == %d; lost_chhild == %d\n",
					tp,
					chlds[tp].pid,
					chlds[tp].status,
					lost_child);
			}
		} while (0);
	}
	
	return 0;
}
