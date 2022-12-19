
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

int main()
{
    const int NX = 400, NY = 400;
    const float DH = 2.0;

    const float vp1 = 3369.0, vs1 = 1643.0, rho1 = 2000.0, Qp1 = 8, Qs1 = 8;
    const float epsilon1 = 0.26, delta1 = -0.05, theta1 = 30.0;

    /*ouput file */
    const char MFILE[50] = "../model/TTI_Thomsen_1";

    float rhov, vp, vs, qp, qs, epsilon, delta, theta;
    int i, j;
    FILE *fp_vs, *fp_vp, *fp_rho, *fp_qp, *fp_qs, *fp_delta, *fp_epsilon, *fp_theta;
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

    sprintf(filename, "%s.epsilon", MFILE);
    fp_epsilon = fopen(filename, "wb");
    if (fp_epsilon == NULL) {
        printf(" Could not open model file for Epsilon values!\n");
	exit(EXIT_FAILURE);
    }

    sprintf(filename, "%s.delta", MFILE);
    fp_delta = fopen(filename, "wb");
    if (fp_delta == NULL) {
        printf(" Could not open model file for Delta values!\n");
	exit(EXIT_FAILURE);
    }

    sprintf(filename, "%s.theta", MFILE);
    fp_theta = fopen(filename, "wb");
    if (fp_theta == NULL) {
        printf(" Could not open model file for Theta values!\n");
	exit(EXIT_FAILURE);
    }

    for (i = 1; i <= NX; i++) {
        for (j = 1; j <= NY; j++) {

            vp = vp1;
            vs = vs1;
            rhov = rho1;
            qp = Qp1;
            qs = Qs1;
            epsilon = epsilon1;
            delta = delta1;
            theta = theta1;

            fwrite(&vp, sizeof(float), 1, fp_vp);
            fwrite(&vs, sizeof(float), 1, fp_vs);
            fwrite(&rhov, sizeof(float), 1, fp_rho);
            fwrite(&qp, sizeof(float), 1, fp_qp);
            fwrite(&qs, sizeof(float), 1, fp_qs);
            fwrite(&epsilon, sizeof(float), 1, fp_epsilon);
            fwrite(&delta, sizeof(float), 1, fp_delta);
            fwrite(&theta, sizeof(float), 1, fp_theta);
        }
    }

    fclose(fp_vp);
    fclose(fp_vs);
    fclose(fp_rho);
    fclose(fp_qp);
    fclose(fp_qs);
    fclose(fp_epsilon);
    fclose(fp_delta);
    fclose(fp_theta);

    sprintf(filename, "%s.vp", MFILE);
    printf("Use for instance \n");
    printf("  ximage n1=%d d1=%f d2=%f < %s label1=\"Y (m)\" label2=\"X (m)\" title=\"%s\"\n", 
	   NY, DH, DH, filename, filename);
    printf("to visualize Vp model.\n");

    return EXIT_SUCCESS;
}
