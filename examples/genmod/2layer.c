
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

int main()
{
    const int NX = 300, NY = 150;
    const float DH = 2.0;

    /* parameters for layer 1 */
    const float vp1 = 2000.0, vs1 = 1000.0, rho1 = 2000.0, Qp1 = 200, Qs1 = 200;

    /* parameters for layer 2  */
    const float vp2 = 3000.0, vs2 = 1500.0, rho2 = 2600.0, Qp2 = 500, Qs2 = 500;

    /* depth of layer 2 (in meter) */
    const float d = 150.0;

    /*ouput file */
    const char MFILE[30] = "../model/2layer";

    float rhov, vp, vs, y, qp, qs;
    int i, j;
    FILE *fp_vs, *fp_vp, *fp_rho, *fp_qp, *fp_qs;
    char filename[75];

    sprintf(filename, "%s.vp", MFILE);
    fp_vp = fopen(filename, "wb");
    if (fp_vp == NULL) {
        printf(" Could not open model file for P-velocities!\n");
	exit(EXIT_FAILURE);
    }

    sprintf(filename, "%s.vs", MFILE);
    fp_vs = fopen(filename, "wb");
    if (fp_vs == NULL) {
        printf(" Could not open model file for S-velocities!\n");
	exit(EXIT_FAILURE);
    }

    sprintf(filename, "%s.rho", MFILE);
    fp_rho = fopen(filename, "wb");
    if (fp_rho == NULL) {
        printf(" Could not open model file for densities!\n");
	exit(EXIT_FAILURE);
    }

    sprintf(filename, "%s.qp", MFILE);
    fp_qp = fopen(filename, "wb");
    if (fp_qp == NULL) {
        printf(" Could not open model file for Qp values!\n");
	exit(EXIT_FAILURE);
    }

    sprintf(filename, "%s.qs", MFILE);
    fp_qs = fopen(filename, "wb");
    if (fp_qs == NULL) {
        printf(" Could not open model file for Qs values!\n");
	exit(EXIT_FAILURE);
    }

    for (i = 1; i <= NX; i++)
        for (j = 1; j <= NY; j++) {

            y = (float)j * DH;

            vp = vp1;
            vs = vs1;
            rhov = rho1;
            qp = Qp1;
            qs = Qs1;

            if (y > d) {
                vp = vp2;
                vs = vs2;
                rhov = rho2;
                qp = Qp2;
                qs = Qs2;
            }

            fwrite(&vp, sizeof(float), 1, fp_vp);
            fwrite(&vs, sizeof(float), 1, fp_vs);
            fwrite(&rhov, sizeof(float), 1, fp_rho);
            fwrite(&qp, sizeof(float), 1, fp_qp);
            fwrite(&qs, sizeof(float), 1, fp_qs);
        }

    fclose(fp_vp);
    fclose(fp_vs);
    fclose(fp_rho);
    fclose(fp_qp);
    fclose(fp_qs);

    sprintf(filename, "%s.vp", MFILE);
    printf("Use for instance \n");
    printf("  ximage n1=%d d1=%f d2=%f < %s label1=\"Y (m)\" label2=\"X (m)\" title=\"%s\"\n", 
	   NY, DH, DH, filename, filename);
    printf("to visualize Vp model.\n");

    return EXIT_SUCCESS;
}
