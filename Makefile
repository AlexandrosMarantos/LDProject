
# Makefile for PatternMatch

CC = gcc

CFLAGS = -O3 -Wall -g -msse4.2 #-D_IPOPCNT 

LIBRARIES = -lm 

EXECNAME = PatternMatch

OBJS = PM.o PM_SupportFunctions.o PM_CommandLine.o PM_ExecutionDetails.o PM_Input.o PM_PatternPool.o PM_Window.o PM_VCFParser.o PM_LinkageDisequilibrium.o PM_LDClass.o

all: $(EXECNAME)

PatternMatch : $(OBJS)
	$(CC) $(CFLAGS) -o $(EXECNAME) $(OBJS) $(LIBRARIES)

PM.o: PM.c 
	$(CC) $(CFLAGS) -c PM.c

PM_SupportFunctions.o: PM_SupportFunctions.c 
	$(CC) $(CFLAGS) -c PM_SupportFunctions.c

PM_CommandLine.o: PM_CommandLine.c 
	$(CC) $(CFLAGS) -c PM_CommandLine.c

PM_ExecutionDetails.o: PM_ExecutionDetails.c 
	$(CC) $(CFLAGS) -c PM_ExecutionDetails.c

PM_Input.o: PM_Input.c 
	$(CC) $(CFLAGS) -c PM_Input.c

PM_PatternPool.o: PM_PatternPool.c 
	$(CC) $(CFLAGS) -c PM_PatternPool.c

PM_Window.o: PM_Window.c 
	$(CC) $(CFLAGS) -c PM_Window.c

PM_VCFParser.o: PM_VCFParser.c 
	$(CC) $(CFLAGS) -c PM_VCFParser.c

PM_LinkageDisequilibrium.o: PM_LinkageDisequilibrium.c 
	$(CC) $(CFLAGS) -c PM_LinkageDisequilibrium.c

PM_LDClass.o: PM_LDClass.c 
	$(CC) $(CFLAGS) -c PM_LDClass.c

clean:
	rm $(EXECNAME)
	rm $(OBJS)
