// pocketflow.cpp
// Author- Pawan Sutar
// Simulation of fluid interaction with WR Wall-Pockets

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define Lx 1.0
#define nx 40
#define Ly 1.5
#define ny 60
#define cfl 0.05
#define U 10.0

int main(void) {
	double u[nx][ny+1], un[nx][ny+1], uc[nx][ny+1], un2[nx][ny+1];
	double v[nx+1][ny], vn[nx+1][ny], vc[nx+1][ny], vn2[nx+1][ny];
	double p[nx + 1][ny + 1], pn[nx + 1][ny + 1], pc[nx + 1][ny + 1], pn2[nx + 1][ny + 1];
	double m[nx + 1][ny + 1];

	int i, j, step;
	double dx, dy, dt, tau, delta, error, Re;
	const int maxstep = 1000000;
	double residual[maxstep], umomresidual[maxstep], vmomresidual[maxstep];
	step = 1;
	dx = Lx / nx;
	dy = Ly / ny;
	dt = cfl * dx / U;
	delta = 1.0;
	error = 1.0;
	Re = 3200.0;
	clock_t start, end;
	double cpu_time_used;
	//double WU = 0.0;

	// Initializing u
	for (i = 0; i <= nx - 1; i++) {
		for (j = 0; j <= ny; j++) {
			u[i][j] = 0.0;
			u[i][ny] = U;
			u[i][ny - 1] = U;
		}
	}

	// Initializing v
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny - 1; j++) {
			v[i][j] = 0.0;
		}
	}

	// Initializing p
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++) {
			p[i][j] = 1.0;
		}
	}

	start = clock();
	while (error < 0.000001) {
		// Solve u-mom
		for (i = 1; i <= nx - 2; i++) {
			for (j = 1; j <= ny - 1; j++) {
				un[i][j] = u[i][j] - dt * ((u[i + 1][j] * u[i + 1][j] - u[i - 1][j] * u[i - 1][j]) / 2.0 / dx
					+ 0.25 * ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j]) - (u[i][j] + u[i][j - 1]) * (v[i + 1][j - 1] + v[i][j - 1])) / dy)
					- dt / dx * (p[i + 1][j] - p[i][j])
					+ dt * 1.0 / Re * ((u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / dx / dx + (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / dy / dy);
			}
		}

		for (j = 1; j = ny - 1; j++) {
			un[0][j] = 0.0;
			un[nx - 1][j] = 0.0;
		}

		for (i = 0; i = nx - 1; i++) {
			un[i][0] = -un[i][1];
			un[i][ny] = 2.0 - un[i][ny - 1];
		}

		// Solve v-mom
		for (i = 0; i = nx - 1; i++) {
			for (j = 0; j = ny - 2; j++) {
				vn[i][j] = v[i][j] - dt * (0.25 * ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j]) - (u[i - 1][j] + u[i - 1][j + 1]) * (v[i][j] + v[i - 1][j])) / dx
					+ (v[i][j + 1] * v[i][j + 1] - v[i][j - 1] * v[i][j - 1]) / 2.0 / dy)
					- dt / dy * (p[i][j + 1] - p[i][j])
					+ dt * 1.0 / Re * ((v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) / dx / dx + (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / dy / dy);
			}
		}

		for (j = 1; j = ny - 2; j++) {
			vn[0][j] = -vn[1][j];
			vn[nx][j] = -vn[nx - 1][j];
		}

		for (i = 0; i <= nx; i++) {
			vn[i][0] = 0.0;
			vn[i][ny - 1] = 0.0;
		}

		// Solve Continuity Equation
		for (i = 1; i <= nx - 1; i++) {
			for (j = 1; j <= ny - 1; j++) {
				pn[i][j] = p[i][j] - dt * delta * ((un[i][j] - un[i - 1][j]) / dx + (vn[i][j] - vn[i][j - 1]) / dy);
			}
		}

		for (i = 1; i <= nx - 1; i++) {
			pn[i][0] = pn[i][1];
			pn[i][ny] = pn[i][ny - 1];
		}

		for (j = 1; j <= ny; j++) {
			pn[0][j] = pn[1][j];
			pn[nx][j] = pn[nx - 1][j];
		}

		// Display Error
		error = 0.0;
		for (i = 0; i <= nx - 1; i++) {
			for (j = 1; j <= ny - 1; j++) {
				m[i][j] = ((un[i][j] - un[i - 1][j]) / dx + (vn[i][j] - vn[i][j - 1]) / dy);
				error = error + fabs(m[i][j]);
			}
		}

		// residual[step] = log10(error)

		if (step % 10000 == 1) {
			printf("Error is %5.101f for step %d\n", error, step);
		}

		// Iterating u
		for (i = 0; i <= nx - 1; i++) {
			for (j = 0; j <= ny; j++) {
				u[i][j] = un[i][j];
			}
		}

		// Iterating v
		for (i = 0; i <= nx; i++) {
			for (j = 0; j <= ny - 1; j++) {
				v[i][j] = vn[i][j];
			}
		}

		// Iterating p
		for (i = 0; i <= nx; i++) {
			for (j = 0; j <= ny; j++) {
				p[i][j] = pn[i][j];
			}
		}

		step = step + 1;
		//WU = WU + 1;

	}

	end = clock();
	cpu_time_used = ((double)(end - start));

	for (i = 0; i <= nx - 1; i++) {
		for (j = 0; j <= ny - 1; j++) {
			uc[i][j] = 0.5 * (u[i][j] + u[i][j + 1]);
			vc[i][j] = 0.5 * (v[i][j] + v[i + 1][j]);
			pc[i][j] = 0.25 * (p[i][j] + p[i + 1][j] + p[i][j + 1] + p[i + 1][j + 1]);
		}
	}

	printf("CPU Time = %lfs\n", cpu_time_used);
	//printf("Work Units = %lf\n", WU);

	// OUTPUT Data

	FILE* fout2, * fout3, * fout4, * fout5, * fout6, * fout7;
	fout2 = fopen("UVP.plt", "w+t");
	fout3 = fopen("Central_U.plt", "w+t");
	fout7 = fopen("Central_V.plt", "w+t");
	fout4 = fopen("Residual.plt", "w+t");
	fout5 = fopen("ResidualUMom.plt", "w+t");
	fout6 = fopen("ResidualVMom.plt", "w+t");

	if (fout2 == NULL) {
		printf("\nError when opening file \n");
		fclose(fout2);
	}
	else {
		fprintf( fout2, "VARIABLES =\"X\", \"Y\", \"U\", \"V\", \"P\"\n");
		fprintf(fout2, "ZONE F=POINT\n");
		fprintf(fout2, "I=%d, J=%d\n", nx, ny);

		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				double xpos, ypos;
				xpos = i * dx;
				ypos = j * dy;
				fprintf(fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, uc[i][j], vc[i][j], pc[i][j]);
			}
		}
	}

	fclose(fout2);

	// Central --U
	fprintf(fout3, "VARIABLES =\"U\", \"Y\"\n");
	fprintf(fout3, "ZONE F=POINT\n");
	fprintf(fout3, "I=%d\n", ny);

	for (j = 0; j < ny; j++) {
		double ypos;
		ypos = (double) j * dy;
		fprintf(fout3, "%5.8lf\t%5.8lf\n", (uc[(nx - 1) / 2][j]), ypos);
	}

	// Central --V
	fprintf(fout7, "VARIABLES =\"V\", \"X\"\n");
	fprintf(fout7, "ZONE F=POINT\n");
	fprintf(fout7, "I=%d\n", nx);

	for (i = 0; i < nx; i++) {
		double xpos;
		xpos = (double) i * dx;
		fprintf(fout7, "%5.8lf\t%5.8lf\n", (vc[i][(ny - 1) / 2]), xpos);
	}



}