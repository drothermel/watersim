#include "Water.h"
#include "MD.h"

int main(int argc, char *argv[])
{
	MD* moldyn = new MD(19.0, 225, 3, 100); //(boxL, nmols, T, steps)
	printf("Made it to the end\n");
	fflush(stdout);
}
