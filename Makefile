CC 		= gcc
CFLAGS  = -g
LIBS    = -lm

SRCS = RACIPE.c pcg_basic.c RACIPELIB.c
OBJS = $(SRCS:.c=.o)

MAIN = RACIPE

.PHONY: depend clean

all: $(MAIN)

$(MAIN):$(OBJS)
		$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^