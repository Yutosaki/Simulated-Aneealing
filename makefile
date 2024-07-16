TARGET = tsp

SRCS = $(TARGET).c cpu_time.c

OBJS = $(SRCS:.c=.o)

CC= gcc
CFLAGS= -Wall -O2 

all: $(TARGET)

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).o -lm

$(TARGET).o: $(TARGET).c cpu_time.c
	$(CC) $(CFLAGS) -c $(TARGET).c

clean:
	rm -f $(OBJS)

fclean: clean
	rm -f $(TARGET)

re: fclean all

.PHONY: all clean fclean re
