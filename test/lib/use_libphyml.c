/* A program that uses the installed library (phyml.h and libphyml.a) the
   way an application would: built and run by `make -C test check-install`
   against `make install`'s files, with nothing from the source tree. */

#include <stdio.h>
#include <string.h>

#include <phyml.h>

static int failures = 0;

static void expect(int cond, const char *what)
{
  printf("%s  %s\n", cond ? "ok  " : "FAIL", what);
  if (!cond) failures++;
}

int main(int argc, char **argv)
{
  char  tree_file[4096];
  FILE *fp;
  int   rc;

  if (argc != 3)
  {
    fprintf(stderr, "usage: %s alignment.phy malformed.phy\n", argv[0]);
    return 2;
  }

  // An error comes back as a code, with its text, and the process goes on
  {
    char *args[] = {"phyml", "-i", "/no/such/alignment.phy", NULL};
    rc           = PhyML_Run(3, args);
    expect(rc == PHYML_ERR_USAGE, "missing file: PHYML_ERR_USAGE");
    expect(strstr(PhyML_Error_Message(), "does not exist") != NULL,
           "missing file: error text");
  }
  {
    char *args[] = {"phyml", "-i", argv[2], "-b", "0", NULL};
    rc           = PhyML_Run(5, args);
    expect(rc == PHYML_ERR_INPUT, "malformed alignment: PHYML_ERR_INPUT");
    expect(!strcmp(PhyML_Error_Name(rc), "PHYML_ERR_INPUT"), "error name");
  }

  // Then a complete analysis
  {
    char *args[] = {"phyml", "-i", argv[1], "-m", "HKY85", "-b", "0",
                    "--r_seed", "1", NULL};
    rc           = PhyML_Run(9, args);
    expect(rc == PHYML_OK, "analysis: PHYML_OK");
  }
  snprintf(tree_file, sizeof(tree_file), "%s_phyml_tree.txt", argv[1]);
  fp = fopen(tree_file, "r");
  expect(fp != NULL && fgetc(fp) == '(', "analysis: tree file written");
  if (fp) fclose(fp);

  printf("%s\n", failures ? "FAILED" : "OK");
  return failures ? 1 : 0;
}
