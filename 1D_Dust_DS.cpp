
#include "pch.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Data.h"



double W(double r, double h)
{
	double k, c, tmp1;
	k = (double) fabs(r) / h;
	c = (double) 2 / 3;
	if (k <= 1.0) { tmp1 = (double) 1.0 - 1.5*k*k + 0.75 * k*k*k; }
	if ((k >= 1.0) && (k <= 2.0)) { tmp1 = (double) 0.25 * (2 - k)*(2 - k)*(2 - k); }
	if (k >= 2.0) { tmp1 = 0.0; }

	return  c / h * tmp1;
}

double dW(double r, double h)
{
	double k, c, tmp1;
	k = (double)r / h;
	c = (double)2 / 3;
	if (k <= -2.0) { tmp1 = 0.0; }
	if ((k >= -2.0) && (k <= -1.0)) { tmp1 = (double) 0.75 * (2.0 + k)*(2.0 + k); }
	if ((k >= -1.0) && (k <= 0)) { tmp1 = (double)-3.0*k - 2.25 * k*k; }
	if ((k >= 0) && (k <= 1.0)) { tmp1 = (double)-3.0*k + 2.25 * k*k; }
	if ((k >= 1.0) && (k <= 2.0)) { tmp1 = (double)-0.75 * (2.0 - k)*(2.0 - k); }
	if (k >= 2.0) { tmp1 = 0.0; }

	return  c / (h*h) * tmp1;
}


double Eq_State(double rho, double e)
{
	//	return B*(pow(rho,Gam_liq)-1) + 1;
	if (Type_of_state == 0) { return rho * e* (Gam_g - 1); }
	if (Type_of_state == 1) { return c_gas * c_gas*rho; }
}

double Eq_State1(double p, double e)
{
	//	return pow((p-1)/B+1,1.0/Gam_liq);

	if (Type_of_state == 0) { return p / (e*(Gam_g - 1)); }
	if (Type_of_state == 1) { return p / (c_gas*c_gas); }
}

double Eq_t_stop(int l, double p, double rho, double gamma)
{
	if (Variant_type == 0) { return R_dust[l] / (sqrt(gamma*p / rho)*rho); }
	if (Variant_type == 1) { return t_stop_init[l]; }
}

double f_init_rho(double x1, double x, double a_cos, double b_sin, double masgas)
{
	// double rho_init = 1.0;
	return rho_init * masgas *(x - x1) + delta * (a_cos*sin(2.0*PI*x) - b_sin * cos(2.0*PI*x) - a_cos * sin(2.0*PI*x1) + b_sin * cos(2.0*PI*x1)) / (2.0*PI) - rho_init * masgas / Num_gas_particle_0;
}

double eq_dihot(double a, double b, double e, double a_cos, double b_sin, double masgas)
{
	double x1, tmp;

	x1 = a;
	tmp = a;
	while ((abs(f_init_rho(x1, a, a_cos, b_sin, masgas) - f_init_rho(x1, b, a_cos, b_sin, masgas)) > e) && (abs(a - b) > e))
	{
		tmp = (b + a) / 2.0;
		if (f_init_rho(x1, a, a_cos, b_sin, masgas)*f_init_rho(x1, tmp, a_cos, b_sin, masgas) < 0.0) { b = tmp; }
		else { a = tmp; }
	}
	return tmp;
}
void Data_out(int num)
{
	FILE  *out_file_gas, *out_file_dust;
	char out_name[25];
	int i, j, l;


	if (num < 10000000) { sprintf(out_name, "Data/G%d.dat", num); };
	if (num < 1000000) { sprintf(out_name, "Data/G0%d.dat", num); };
	if (num < 100000) { sprintf(out_name, "Data/G00%d.dat", num); };
	if (num < 10000) { sprintf(out_name, "Data/G000%d.dat", num); };
	if (num < 1000) { sprintf(out_name, "Data/G0000%d.dat", num); };
	if (num < 100) { sprintf(out_name, "Data/G00000%d.dat", num); };
	if (num < 10) { sprintf(out_name, "Data/G000000%d.dat", num); };

	out_file_gas = fopen(out_name, "wt");
	fprintf(out_file_gas, "t=%5.3f \n", Tm);
	fprintf(out_file_gas, "tau=%10.8lf \t h=%10.8lf \t h_cell=%10.8lf \n", tau, h0, cell_width);
	fprintf(out_file_gas, "x \t mas \t rho \t p \t v \t acc \t e \t Ind \n");

	for (i = 0; i <= Num_gas_particle; i++)
	{
		if ((x_gas[i] >= 0.0) && (x_gas[i] <= 1.0))
		{
			fprintf(out_file_gas, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %d \n",
				x_gas[i], mas_gas[i], rho_gas[i], p_gas[i], vel_gas[i], acc_gas[i], e_gas[i], ind_gas[i]);
		}

	}

	fclose(out_file_gas);

	for (l = 0; l < Num_dust_sort; l++)
	{
		if (num < 10000000) { sprintf(out_name, "Data/D-%d-%d.dat", l, num); };
		if (num < 1000000) { sprintf(out_name, "Data/D-%d-0%d.dat", l, num); };
		if (num < 100000) { sprintf(out_name, "Data/D-%d-00%d.dat", l, num); };
		if (num < 10000) { sprintf(out_name, "Data/D-%d-000%d.dat", l, num); };
		if (num < 1000) { sprintf(out_name, "Data/D-%d-0000%d.dat", l, num); };
		if (num < 100) { sprintf(out_name, "Data/D-%d-00000%d.dat", l, num); };
		if (num < 10) { sprintf(out_name, "Data/D-%d-000000%d.dat", l, num); };

		out_file_dust = fopen(out_name, "wt");
		fprintf(out_file_dust, "t=%5.3f \n", Tm);
		fprintf(out_file_dust, "x \t mas \t rho \t v \t acc \t t_stop \t Ind \n");

		for (j = 0; j <= Num_dust_particle; j++)
		{
			if ((x_dust[l][j] >= 0.0) && (x_dust[l][j] <= 1.0))
			{
				fprintf(out_file_dust, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %d \n",
					x_dust[l][j], mas_dust[l][j], rho_dust[l][j], vel_dust[l][j], acc_dust[l][j], t_stop[l][j], ind_dust[l][j]);
			}
		}

		fclose(out_file_dust);
	}

}

int main()
{
	int  i, j, i1, j1, I1l, I1r, l;
	char s[128];
	double K, h_left, h_right, e_left, e_right, Coeff_h_cell;
	double tmp1, tmp2, tmp3;

	FILE *in_file, *out_particle_file;
	char out_name[8];
	int out_num;
	int dh = 50;
	int cell_num;
	int Num_out_particle;
	bool Flag_num_out_particle;

	out_particle_file = nullptr;

	printf("Dust \n");


	// ���� ������

	in_file = fopen("Init.txt", "rt");

	fgets(s, 128, in_file);
	fgets(s, 128, in_file);	sscanf(s, "%d", &Variant_type);
	fgets(s, 128, in_file);	sscanf(s, "%d", &Flag_of_method);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &Rm);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &c0);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &c_gas);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &Gam_liq);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &Gam_g);
	fgets(s, 128, in_file); sscanf(s, "%lf", &p_left);
	fgets(s, 128, in_file); sscanf(s, "%lf", &p_right);
	fgets(s, 128, in_file);	sscanf(s, "%d", &Num_gas_particle_0);
	fgets(s, 128, in_file); sscanf(s, "%lf", &tau);
	fgets(s, 128, in_file); sscanf(s, "%lf", &h0);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &T_end);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &T_out);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &Coeff_h_cell);
	fgets(s, 128, in_file); sscanf(s, "%lf", &alpha);
	fgets(s, 128, in_file); sscanf(s, "%lf", &beta);
	fgets(s, 128, in_file);	sscanf(s, "%d", &Num_dust_sort);
	for (l = 0; l < Num_dust_sort; l++)
	{
		fgets(s, 128, in_file); sscanf(s, "%lf", &R_dust[l]);
		fgets(s, 128, in_file); sscanf(s, "%lf", &mas_dustgas[l]);
	}
	fclose(in_file);


	dh = (int)(2.0*h0 / (1.0 / Num_gas_particle_0)) + 5;
	mas_gas_0 = 1.0 / Num_gas_particle_0;
	rho_init = 1.0;

	eps = 0.01;

	if (Variant_type == 0)
	{
		Type_of_state = 0;
		e_left = 2.5;
		e_right = 2.0;
		rho_left = Eq_State1(p_left, e_left);
		rho_right = Eq_State1(p_right, e_right);

		Xmin = -0.2*Rm;
		Xmax = 1.2 * Rm;
		Num_gas_particle = (int)((Xmax - Xmin)*0.5 * rho_left / mas_gas_0 + (Xmax - Xmin)*0.5 * rho_right / mas_gas_0) + 1;

		Num_dust_particle = Num_gas_particle;

		h_left = mas_gas_0 / rho_left;
		h_right = mas_gas_0 / rho_right;

		x_gas = new double[Num_gas_particle + 1];
		h_gas = new double[Num_gas_particle + 1];
		p_gas = new double[Num_gas_particle + 1];
		mas_gas = new double[Num_gas_particle + 1];
		rho_gas = new double[Num_gas_particle + 1];
		vel_gas = new double[Num_gas_particle + 1];
		acc_gas = new double[Num_gas_particle + 1];
		e_gas = new double[Num_gas_particle + 1];
		ind_gas = new int[Num_gas_particle + 1];


		for (l = 0; l < Num_dust_sort; l++)
		{
			x_dust[l] = new double[Num_gas_particle + 1];
			h_dust[l] = new double[Num_gas_particle + 1];
			mas_dust[l] = new double[Num_gas_particle + 1];
			rho_dust[l] = new double[Num_gas_particle + 1];
			vel_dust[l] = new double[Num_gas_particle + 1];
			acc_dust[l] = new double[Num_gas_particle + 1];
			ind_dust[l] = new int[Num_gas_particle + 1];
			t_stop[l] = new double[Num_gas_particle + 1];
		}


		Im1_gas = dh;
		Im2_gas = Num_gas_particle - dh;

		x_gas[0] = Xmin;
		vel_gas[0] = 0.0;
		acc_gas[0] = 0.0;
		h_gas[0] = h0;
		p_gas[0] = p_left;
		rho_gas[0] = rho_left;
		e_gas[0] = e_left;
		mas_gas[0] = mas_gas_0;
		ind_gas[0] = 10;

		for (i = 1; i <= Num_gas_particle; i++)
		{
			vel_gas[i] = 0.0;
			acc_gas[i] = 0.0;
			h_gas[i] = h0;
			if (i < (int)((Xmax - Xmin)*0.5 * rho_left / mas_gas_0))
			{
				x_gas[i] = x_gas[i - 1] + h_left;
				p_gas[i] = p_left;
				rho_gas[i] = rho_left;
				e_gas[i] = e_left;
				mas_gas[i] = mas_gas_0;
			}
			else
			{
				x_gas[i] = x_gas[i - 1] + h_right;
				p_gas[i] = p_right;
				rho_gas[i] = rho_right;
				e_gas[i] = e_right;
				mas_gas[i] = mas_gas_0;
			};
			ind_gas[i] = 0;
			if ((i < Im1_gas) || (i > Im2_gas)) { ind_gas[i] = 10; }
		}


		Im1_dust = dh;
		Im2_dust = Num_dust_particle - dh;
		for (l = 0; l < Num_dust_sort; l++)
		{
			rho_dust_left = rho_left * mas_dustgas[l];
			rho_dust_right = rho_right * mas_dustgas[l];
			t_stop[l][0] = Eq_t_stop(l, p_gas[0], rho_gas[0], Gam_g);

			x_dust[l][0] = Xmin;
			vel_dust[l][0] = 0.0;
			acc_dust[l][0] = 0.0;
			h_dust[l][0] = h0;
			rho_dust[l][0] = rho_dust_left;
			mas_dust[l][0] = mas_gas_0 * mas_dustgas[l];
			ind_dust[l][0] = 10;



			for (j = 1; j <= Num_dust_particle; j++)
			{

				vel_dust[l][j] = 0.0;
				acc_dust[l][j] = 0.0;
				h_dust[l][j] = h0;
				if (j < (int)((Xmax - Xmin)*0.5  * rho_left / mas_gas_0))
				{
					x_dust[l][j] = x_dust[l][j - 1] + h_left;
					rho_dust[l][j] = rho_dust_left;
					mas_dust[l][j] = mas_gas_0 * mas_dustgas[l];
				}
				else
				{
					x_dust[l][j] = x_dust[l][j - 1] + h_right;
					rho_dust[l][j] = rho_dust_right;
					mas_dust[l][j] = mas_gas_0 * mas_dustgas[l];
				};
				ind_dust[l][j] = 0;
				t_stop[l][j] = Eq_t_stop(l, p_gas[j], rho_gas[j], Gam_g);
				if ((j < Im1_dust) || (j > Im2_dust)) { ind_dust[l][j] = 10; }
			}
		}

		Number_of_cell = (int)((Xmax - Xmin) / (Coeff_h_cell* h0));
		x_cell = new double[Number_of_cell + 1];
		v_g_average = new double[Number_of_cell + 1];
		rho_g_average = new double[Number_of_cell + 1];
		e_g_average = new double[Number_of_cell + 1];
		g_average_count = new int[Number_of_cell + 1];
		Psi_av_new = new double[Number_of_cell + 1];
		v_av_new = new double[Number_of_cell + 1];
		y_av_new = new double[Number_of_cell + 1];
		y_av = new double[Number_of_cell + 1];
		beta_cell = new double[Number_of_cell + 1];

		for (l = 0; l < Num_dust_sort; l++)
		{
			v_d_average[l] = new double[Number_of_cell + 1];
			rho_d_average[l] = new double[Number_of_cell + 1];
			eps_cell[l] = new double[Number_of_cell + 1];
			d_average_count[l] = new int[Number_of_cell + 1];
			t_stop_average[l] = new double[Number_of_cell + 1];
			x_av_new[l] = new double[Number_of_cell + 1];
			x_av[l] = new double[Number_of_cell + 1];
			u_av_new[l] = new double[Number_of_cell + 1];
			b_cell[l] = new double[Number_of_cell + 1];
		}

		cell_width = (Xmax - Xmin) / Number_of_cell;
		for (i = 0; i <= Number_of_cell; i++)
		{
			x_cell[i] = (double)Xmin + i * cell_width;
			v_g_average[i] = 0.0;
			rho_g_average[i] = 0.0;
			e_g_average[i] = 0.0;
			g_average_count[i] = 0;
			Psi_av_new[i] = 0.0;
			v_av_new[i] = 0.0;
			y_av_new[i] = 0.0;
			y_av[i] = 0.0;
			for (l = 0; l < Num_dust_sort; l++)
			{
				v_d_average[l][i] = 0.0;
				rho_d_average[l][i] = 0.0;
				eps_cell[l][i] = 0.0;
				t_stop_average[l][i] = 0.0;
				d_average_count[l][i] = 0;
				x_av_new[l][i] = 0.0;
				x_av[l][i] = 0.0;
				u_av_new[l][i] = 0.0;
				b_cell[l][i] = 0.0;
			}
		}
	}

	if (Variant_type == 1)
	{
		alpha = 0.0;
		beta = 0.0;
		Type_of_state = 1;

		t_stop_init[0] = 0.1;
		t_stop_init[1] = 0.2;
		t_stop_init[2] = 0.4;
		t_stop_init[3] = 0.1;

		a_rho_gas = 1.0;  b_rho_gas = 0.0;
		a_rho[0] = 0.2813014; b_rho[0] = 0.1508098;
		a_rho[1] = 0.1667321; b_rho[1] = 0.1957177;
		a_rho[2] = 0.0520914; b_rho[2] = 0.1508957;
		a_rho[3] = 0.4309959; b_rho[3] = 0.1901821;

		a_vel_gas = -0.7852741; b_vel_gas = 0.1267991;
		a_vel[0] = -0.7201359; b_vel[0] = -0.2482996;
		a_vel[1] = -0.4672883; b_vel[1] = -0.3976914;
		a_vel[2] = -0.1801365; b_vel[2] = -0.3357016;
		a_vel[3] = -0.1801365; b_vel[3] = -0.219647;


		Xmin = -3.5*Rm;
		Xmax = 4.5 * Rm;


		Num_gas_particle = (int)((Xmax - Xmin)*Num_gas_particle_0) + 1;

		Num_dust_particle = Num_gas_particle;

		x_gas = new double[Num_gas_particle + 1];
		h_gas = new double[Num_gas_particle + 1];
		p_gas = new double[Num_gas_particle + 1];
		mas_gas = new double[Num_gas_particle + 1];
		rho_gas = new double[Num_gas_particle + 1];
		vel_gas = new double[Num_gas_particle + 1];
		acc_gas = new double[Num_gas_particle + 1];
		e_gas = new double[Num_gas_particle + 1];
		ind_gas = new int[Num_gas_particle + 1];


		for (l = 0; l < Num_dust_sort; l++)
		{
			x_dust[l] = new double[Num_gas_particle + 1];
			h_dust[l] = new double[Num_gas_particle + 1];
			mas_dust[l] = new double[Num_gas_particle + 1];
			rho_dust[l] = new double[Num_gas_particle + 1];
			vel_dust[l] = new double[Num_gas_particle + 1];
			acc_dust[l] = new double[Num_gas_particle + 1];
			ind_dust[l] = new int[Num_gas_particle + 1];
			t_stop[l] = new double[Num_gas_particle + 1];
		}


		Im1_gas = dh;
		Im2_gas = Num_gas_particle - dh;

		delta = 0.0001;

		x_gas[0] = Xmin;
		acc_gas[0] = 0.0;
		h_gas[0] = h0;
		rho_gas[0] = rho_init + delta * (a_rho_gas * cos(2.0*PI*x_gas[0]) + b_rho_gas * sin(2.0*PI*x_gas[0]));
		vel_gas[0] = delta * (a_vel_gas  * cos(2.0*PI*x_gas[0]) + b_vel_gas * sin(2.0*PI*x_gas[0]));
		e_gas[0] = 1.0;
		mas_gas[0] = mas_gas_0;
		p_gas[0] = Eq_State(rho_gas[0], e_gas[0]);
		ind_gas[0] = 10;

		Flag_num_out_particle = true;

		for (i = 1; i <= Num_gas_particle; i++)
		{
			acc_gas[i] = 0.0;
			h_gas[i] = h0;
			x_gas[i] = eq_dihot(x_gas[i - 1], x_gas[i - 1] + 2.0 / Num_gas_particle_0, 1.0e-15, a_rho_gas, b_rho_gas, 1.0);
			rho_gas[i] = rho_init + delta * (a_rho_gas * cos(2.0*PI*x_gas[i]) + b_rho_gas * sin(2.0*PI*x_gas[i]));
			vel_gas[i] = delta * (a_vel_gas  * cos(2.0*PI*x_gas[i]) + b_vel_gas * sin(2.0*PI*x_gas[i]));
			e_gas[i] = 1.0;
			mas_gas[i] = mas_gas_0;
			p_gas[i] = Eq_State(rho_gas[i], e_gas[i]);
			ind_gas[i] = 0;
			if ((i < Im1_gas) || (i > Im2_gas)) { ind_gas[i] = 10; }
			if (Flag_num_out_particle) { Num_out_particle = i; };
			if (x_gas[i] >= 0.0) { Flag_num_out_particle = false; };
		}


		Im1_dust = dh;
		Im2_dust = Num_dust_particle - dh;
		for (l = 0; l < Num_dust_sort; l++)
		{
			t_stop[l][0] = Eq_t_stop(l, p_gas[0], rho_gas[0], Gam_g);
			x_dust[l][0] = Xmin;
			acc_dust[l][0] = 0.0;
			h_dust[l][0] = h0;
			rho_dust[l][0] = rho_init * mas_dustgas[l] + delta * (a_rho[l] * cos(2.0*PI*x_dust[l][0]) + b_rho[l] * sin(2.0*PI*x_dust[l][0]));
			vel_dust[l][0] = delta * (a_vel[l] * cos(2.0*PI*x_dust[l][0]) + b_vel[l] * sin(2.0*PI*x_dust[l][0]));
			mas_dust[l][0] = mas_gas_0 * mas_dustgas[l];
			ind_dust[l][0] = 10;



			for (j = 1; j <= Num_dust_particle; j++)
			{

				acc_dust[l][j] = 0.0;
				h_dust[l][j] = h0;
				x_dust[l][j] = eq_dihot(x_dust[l][j - 1], x_dust[l][j - 1] + 2.0 / Num_gas_particle_0, 1.0e-15, a_rho[l], b_rho[l], mas_dustgas[l]);
				rho_dust[l][j] = rho_init * mas_dustgas[l] + delta * (a_rho[l] * cos(2.0*PI*x_dust[l][j]) + b_rho[l] * sin(2.0*PI*x_dust[l][j]));
				vel_dust[l][j] = delta * (a_vel[l] * cos(2.0*PI*x_dust[l][j]) + b_vel[l] * sin(2.0*PI*x_dust[l][j]));
				mas_dust[l][j] = mas_gas_0 * mas_dustgas[l];
				ind_dust[l][j] = 0;
				t_stop[l][j] = Eq_t_stop(l, p_gas[j], rho_gas[j], Gam_g);
				if ((j < Im1_dust) || (j > Im2_dust)) { ind_dust[l][j] = 10; }
			}
		}

		Number_of_cell = (int)((Xmax - Xmin) / (Coeff_h_cell* h0));
		x_cell = new double[Number_of_cell + 1];
		v_g_average = new double[Number_of_cell + 1];
		rho_g_average = new double[Number_of_cell + 1];
		e_g_average = new double[Number_of_cell + 1];
		g_average_count = new int[Number_of_cell + 1];
		Psi_av_new = new double[Number_of_cell + 1];
		v_av_new = new double[Number_of_cell + 1];
		y_av_new = new double[Number_of_cell + 1];
		y_av = new double[Number_of_cell + 1];
		beta_cell = new double[Number_of_cell + 1];

		for (l = 0; l < Num_dust_sort; l++)
		{
			v_d_average[l] = new double[Number_of_cell + 1];
			rho_d_average[l] = new double[Number_of_cell + 1];
			eps_cell[l] = new double[Number_of_cell + 1];
			d_average_count[l] = new int[Number_of_cell + 1];
			t_stop_average[l] = new double[Number_of_cell + 1];
			x_av_new[l] = new double[Number_of_cell + 1];
			x_av[l] = new double[Number_of_cell + 1];
			u_av_new[l] = new double[Number_of_cell + 1];
			b_cell[l] = new double[Number_of_cell + 1];
		}

		cell_width = (Xmax - Xmin) / Number_of_cell;
		for (i = 0; i <= Number_of_cell; i++)
		{
			x_cell[i] = (double)Xmin + i * cell_width;
			v_g_average[i] = 0.0;
			rho_g_average[i] = 0.0;
			e_g_average[i] = 0.0;
			g_average_count[i] = 0;
			Psi_av_new[i] = 0.0;
			v_av_new[i] = 0.0;
			y_av_new[i] = 0.0;
			y_av[i] = 0.0;
			for (l = 0; l < Num_dust_sort; l++)
			{
				v_d_average[l][i] = 0.0;
				rho_d_average[l][i] = 0.0;
				eps_cell[l][i] = 0.0;
				t_stop_average[l][i] = 0.0;
				d_average_count[l][i] = 0;
				x_av_new[l][i] = 0.0;
				x_av[l][i] = 0.0;
				u_av_new[l][i] = 0.0;
				b_cell[l][i] = 0.0;
			}
		}



	}


	Tm = 0.0;
	out_num = 0;

	Data_out(out_num);
	out_num = out_num + 1;

	if (Variant_type == 1)
	{
		out_particle_file = fopen("Data/Particle.dat", "wt");
		fprintf(out_particle_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t", Tm, x_gas[Num_out_particle], rho_gas[Num_out_particle] * 0.3333, vel_gas[Num_out_particle]);
		for (l = 0; l < Num_dust_sort; l++)
		{

			fprintf(out_particle_file, "%10.8lf \t %10.8lf \t", rho_dust[l][Num_out_particle], vel_dust[l][Num_out_particle]);
		}
		fprintf(out_particle_file, "\n");
	}

	do
	{

		Tm = Tm + tau;
		printf("Time %5.3lf mks \n", Tm*1e3);

		for (i = 0; i <= Num_gas_particle; i++)
		{
			//		rho_p1[i] = rho_gas[i];
			if (ind_gas[i] == 0)
			{
				tmp1 = 0.0;
				I1l = i - dh;
				I1r = i + dh;
				if (I1l < 0) { I1l = 0; }
				if (I1r > Num_gas_particle) { I1r = Num_gas_particle; }
				for (i1 = I1l; i1 <= I1r; i1++)
				{
					if (ind_gas[i1] % 10 == 0) {
						tmp1 = tmp1 + mas_gas[i1] * W(x_gas[i1] - x_gas[i], h_gas[i]);
					}
				}
				rho_gas[i] = tmp1;
			}
		}

		for (l = 0; l < Num_dust_sort; l++)
		{
			for (j = 0; j <= Num_gas_particle; j++)
			{
				if (ind_dust[l][j] == 0)
				{
					tmp1 = 0.0;
					I1l = j - dh;
					I1r = j + dh;
					if (I1l < 0) { I1l = 0; }
					if (I1r > Num_dust_particle) { I1r = Num_dust_particle; }
					for (i1 = I1l; i1 <= I1r; i1++)
					{
						if (ind_dust[l][i1] % 10 == 0) {
							tmp1 = tmp1 + mas_dust[l][i1] * W(x_dust[l][i1] - x_dust[l][j], h_dust[l][j]);
						}
					}
					rho_dust[l][j] = tmp1;
				}
				t_stop[l][j] = Eq_t_stop(l, p_gas[j], rho_gas[j], Gam_g);
			}
		}



		for (i = 0; i <= Num_gas_particle; i++)
		{
			p_gas[i] = Eq_State(rho_gas[i], e_gas[i]);
			acc_gas[i] = 0.0;
		}


		for (l = 0; l < Num_dust_sort; l++)
		{
			for (j = 0; j <= Num_dust_particle; j++)
			{
				acc_dust[l][j] = 0.0;

			}
		}



		if (Flag_of_method == 2)
		{
			for (i = 0; i <= Number_of_cell; i++)
			{
				v_g_average[i] = 0.0;
				rho_g_average[i] = 0.0;
				e_g_average[i] = 0.0;
				g_average_count[i] = 0;
				Psi_av_new[i] = 0.0;
				v_av_new[i] = 0.0;
				y_av_new[i] = 0.0;
				y_av[i] = 0.0;
			}

			for (l = 0; l < Num_dust_sort; l++)
			{
				for (i = 0; i <= Number_of_cell; i++)
				{
					v_d_average[l][i] = 0.0;
					rho_d_average[l][i] = 0.0;
					eps_cell[l][i] = 0.0;
					d_average_count[l][i] = 0;
					x_av_new[l][i] = 0.0;
					x_av[l][i] = 0.0;
					u_av_new[l][i] = 0.0;
					t_stop_average[l][i] = 0.0;
					b_cell[l][i] = 0.0;
				}
			}

			for (i = 0; i <= Num_gas_particle; i++)
			{
				cell_num = trunc((x_gas[i] - Xmin) / cell_width);
				g_average_count[cell_num] = g_average_count[cell_num] + 1;
				v_g_average[cell_num] = v_g_average[cell_num] + vel_gas[i];
				rho_g_average[cell_num] = rho_g_average[cell_num] + rho_gas[i];
				e_g_average[cell_num] = e_g_average[cell_num] + e_gas[i];
			}

			for (i = 0; i <= Number_of_cell; i++)
			{
				if (g_average_count[i] > 0) { v_g_average[i] = v_g_average[i] / g_average_count[i]; }
				if (g_average_count[i] > 0) { e_g_average[i] = e_g_average[i] / g_average_count[i]; }
			}

			for (l = 0; l < Num_dust_sort; l++)
			{

				for (j = 0; j <= Num_dust_particle; j++)
				{
					cell_num = trunc((x_dust[l][j] - Xmin) / cell_width);
					d_average_count[l][cell_num] = d_average_count[l][cell_num] + 1;
					t_stop_average[l][cell_num] = t_stop_average[l][cell_num] + t_stop[l][j];
					v_d_average[l][cell_num] = v_d_average[l][cell_num] + vel_dust[l][j];
					rho_d_average[l][cell_num] = rho_d_average[l][cell_num] + rho_dust[l][j];
				}
			}
			for (l = 0; l < Num_dust_sort; l++)
			{
				for (i = 0; i <= Number_of_cell; i++)
				{
					if (d_average_count[l][i] > 0) { v_d_average[l][i] = v_d_average[l][i] / d_average_count[l][i]; }
					if (d_average_count[l][i] > 0) { rho_d_average[l][i] = rho_d_average[l][i] / d_average_count[l][i]; }
					if (d_average_count[l][i] > 0) { t_stop_average[l][i] = t_stop_average[l][i] / d_average_count[l][i]; }
					eps_cell[l][i] = d_average_count[l][i] * mas_dustgas[l] / g_average_count[i];
				}
			}


			for (i = 0; i <= Num_gas_particle; i++)
			{
				double tmp1, nu1, Fnu;
				tmp1 = 0.0;
				tmp2 = 0.0;
				j = i;
				if (ind_gas[i] == 0)
				{
					I1l = i - dh;
					I1r = i + dh;
					if (I1l < 0) { I1l = 0; }
					if (I1r > Num_gas_particle) { I1r = Num_gas_particle; }
					for (i1 = I1l; i1 <= I1r; i1++)
					{

						nu1 = (vel_gas[i] - vel_gas[i1])*(x_gas[i] - x_gas[i1]);
						if (nu1 < 0)
						{
							Cnu = 0.5*(sqrt(Gam_g*p_gas[i] / rho_gas[i]) + sqrt(Gam_g*p_gas[i1] / rho_gas[i1]));
							nu1 = h0 * nu1 / ((x_gas[i] - x_gas[i1])*(x_gas[i] - x_gas[i1]) + eps * h0*h0);
							//		        Fnu = 2*(-alfa * Cnu * nu1 + beta * nu1 * nu1) / (p_p[i]+p_p[i1]);
							Fnu = (-alpha * Cnu * nu1 + beta * nu1 * nu1) / (0.5*(rho_gas[i] + rho_gas[i1]));
						}
						else { Fnu = 0; }

						cell_num = trunc((x_gas[i] - Xmin) / cell_width);
						Psi_av_new[cell_num] = Psi_av_new[cell_num] + mas_gas[i1] * (p_gas[i] / (rho_gas[i] * rho_gas[i]) + p_gas[i1] / (rho_gas[i1] * rho_gas[i1]) + Fnu)*dW(x_gas[i1] - x_gas[i], h_gas[i]);
						tmp2 = tmp2 + mas_gas[i] * p_gas[i] / (rho_gas[i] * rho_gas[i]) *(vel_gas[i] - vel_gas[i1])*dW(x_gas[i] - x_gas[i1], h_gas[i]) + mas_gas[i] / 2.0 * Fnu*(vel_gas[i] - vel_gas[i1])*dW(x_gas[i] - x_gas[i1], h_gas[i]);


					}
					e_gas[i] = e_gas[i] + tmp2 * tau;

				}
			}

			for (i = 0; i <= Number_of_cell; i++)
			{
				if (g_average_count[i] > 0) { Psi_av_new[i] = Psi_av_new[i] / g_average_count[i]; }
				y_av[i] = v_g_average[i];
				for (l = 0; l < Num_dust_sort; l++)
				{
					x_av[l][i] = v_g_average[i] - v_d_average[l][i];
					y_av[i] = y_av[i] + eps_cell[l][i] * v_d_average[l][i];
					b_cell[l][i] = (t_stop_average[l][i] + tau) / (eps_cell[l][i] * tau);
				}
			}

			for (i = 0; i <= Number_of_cell; i++)
			{
				beta_cell[i] = 1.0;
				for (l = 0; l < Num_dust_sort; l++)
				{
					beta_cell[i] = beta_cell[i] + 1.0 / b_cell[l][i];
				}
			}

			for (i = 0; i <= Number_of_cell; i++)
			{
				y_av_new[i] = y_av[i] + tau * Psi_av_new[i];
				for (l = 0; l < Num_dust_sort; l++)
				{
					x_av_new[l][i] = -1.0 * b_cell[l][i] * beta_cell[i] * (x_av[l][i] + tau * Psi_av_new[i]) / (b_cell[l][i] * b_cell[l][i]);
					for (j = 0; j < Num_dust_sort; j++)
					{
						x_av_new[l][i] = x_av_new[l][i] + (x_av[j][i] + tau * Psi_av_new[i]) / (b_cell[l][i] * b_cell[j][i]);
					}
					x_av_new[l][i] = x_av_new[l][i] * (-1.0*t_stop_average[l][i]) / (tau * eps_cell[l][i] * beta_cell[i]);
				}
			}

			for (i = 0; i <= Number_of_cell; i++)
			{
				tmp1 = 1.0;
				for (l = 0; l < Num_dust_sort; l++)
				{
					tmp1 = tmp1 + eps_cell[l][i];
				}

				v_g_average[i] = y_av_new[i];

				for (l = 0; l < Num_dust_sort; l++)
				{
					v_g_average[i] = v_g_average[i] + eps_cell[l][i] * x_av_new[l][i];

					v_d_average[l][i] = y_av_new[i];

					for (j = 0; j < Num_dust_sort; j++)
					{
						if (j != l)
						{
							v_d_average[l][i] = v_d_average[l][i] + eps_cell[j][i] * x_av_new[j][i];
						}

					}
					v_d_average[l][i] = v_d_average[l][i] - (tmp1 - eps_cell[l][i]) *  x_av_new[l][i];
					v_d_average[l][i] = v_d_average[l][i] / tmp1;
				}
				v_g_average[i] = v_g_average[i] / tmp1;
			}


			for (l = 0; l < Num_dust_sort; l++)
			{
				for (j = 0; j <= Num_dust_particle; j++)
				{
					cell_num = trunc((x_dust[l][j] - Xmin) / cell_width);
					vel_dust[l][j] = (vel_dust[l][j] / tau + v_g_average[cell_num] / t_stop_average[l][cell_num]) / (1.0 / tau + 1.0 / t_stop_average[l][cell_num]);
				}
			}


			for (i = 0; i <= Num_gas_particle; i++)
			{
				if (ind_gas[i] == 0)
				{
					cell_num = trunc((x_gas[i] - Xmin) / cell_width);
					tmp1 = 0.0;
					tmp2 = 0.0;
					for (l = 0; l < Num_dust_sort; l++)
					{
						tmp1 = tmp1 + eps_cell[l][cell_num] / t_stop_average[l][cell_num];
						tmp2 = tmp2 + eps_cell[l][cell_num] * v_d_average[l][cell_num] / t_stop_average[l][cell_num];
					}

					vel_gas[i] = (vel_gas[i] / tau + tmp2 + Psi_av_new[cell_num]) / (1.0 / tau + tmp1);
				}
			}



			for (i = 0; i <= Num_gas_particle; i++)
			{
				if (ind_gas[i] == 0)
				{
					x_gas[i] = x_gas[i] + vel_gas[i] * tau;
				}
			}

			for (l = 0; l < Num_dust_sort; l++)
			{
				for (j = 0; j <= Num_dust_particle; j++)
				{
					if (ind_dust[l][j] == 0)
					{
						x_dust[l][j] = x_dust[l][j] + vel_dust[l][j] * tau;
					}
				}
			}
		}



		if (Flag_of_method == 1)
		{
			for (i = 0; i <= Num_gas_particle; i++)
			{
				double tmp1, nu1, Fnu;
				tmp1 = 0.0;
				tmp2 = 0.0;
				tmp3 = 0.0;

				if (ind_gas[i] == 0)
				{
					I1l = i - dh;
					I1r = i + dh;
					if (I1l < 0) { I1l = 0; }
					if (I1r > Num_gas_particle) { I1r = Num_gas_particle; }
					for (i1 = I1l; i1 <= I1r; i1++)
					{

						if (ind_gas[i1] % 10 == 0)
						{
							nu1 = (vel_gas[i] - vel_gas[i1])*(x_gas[i] - x_gas[i1]);
							if (nu1 < 0)
							{
								Cnu = 0.5*(sqrt(Gam_g*p_gas[i] / rho_gas[i]) + sqrt(Gam_g*p_gas[i1] / rho_gas[i1]));
								nu1 = h0 * nu1 / ((x_gas[i] - x_gas[i1])*(x_gas[i] - x_gas[i1]) + eps * h0*h0);
								//		        Fnu = 2*(-alfa * Cnu * nu1 + beta * nu1 * nu1) / (p_p[i]+p_p[i1]);
								Fnu = (-alpha * Cnu * nu1 + beta * nu1 * nu1) / (0.5*(rho_gas[i] + rho_gas[i1]));
							}
							else { Fnu = 0; }
							//						K = c0 *abs(vel[i] - vel[j1]) / R_dust;
							//			K = c0 *rho_gas[i] * rho_gas[j1] / R_dust;
							//			tmp2 = mas_gas[j1] * K / (rho_gas[i] * rho_gas[j1]) *(vel_gas[i] - vel_gas[j1])*(x_gas[j1] - x[i])*(x[j1] - x[i]) / ((x[j1] - x[i])*(x[j1] - x[i]) + h[i] * h[i] * 1.0e-6)*
							//				(0.5*W(x[j1] - x[i], h[i]) + 0.5*W(x[j1] - x[i], h[j1]));

							tmp1 = tmp1 + mas_gas[i1] * (p_gas[i] / (rho_gas[i] * rho_gas[i]) + p_gas[i1] / (rho_gas[i1] * rho_gas[i1]) + Fnu)*dW(x_gas[i1] - x_gas[i], h_gas[i]);
							tmp2 = tmp2 + mas_gas[i] * p_gas[i] / (rho_gas[i] * rho_gas[i]) *(vel_gas[i] - vel_gas[i1])*dW(x_gas[i] - x_gas[i1], h_gas[i]) + mas_gas[i] / 2.0 * Fnu*(vel_gas[i] - vel_gas[i1])*dW(x_gas[i] - x_gas[i1], h_gas[i]);

						}
					}
					acc_gas[i] = tmp1;
					e_gas[i] = e_gas[i] + tmp2 * tau;
				}
			}


			for (i = 0; i <= Num_gas_particle; i++)
			{
				if (ind_gas[i] == 0)
				{
					vel_gas[i] = vel_gas[i] + acc_gas[i] * tau;
					x_gas[i] = x_gas[i] + vel_gas[i] * tau;
				}
			}

			for (l = 0; l < Num_dust_sort; l++)
			{
				for (j = 0; j <= Num_dust_particle; j++)
				{
					if (ind_dust[l][j] == 0)
					{
						vel_dust[l][j] = vel_dust[l][j] + acc_dust[l][j] * tau;
						x_dust[l][j] = x_dust[l][j] + vel_dust[l][j] * tau;

					}
				}
			}
		}



		if (Tm >= out_num * T_out)
		{
			Data_out(out_num);
			out_num = out_num + 1;
		}

		if (Variant_type == 1)
		{
			fprintf(out_particle_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t", Tm, x_gas[Num_out_particle], rho_gas[Num_out_particle] * 0.3333, vel_gas[Num_out_particle]);
			for (l = 0; l < Num_dust_sort; l++)
			{

				fprintf(out_particle_file, "%10.8lf \t %10.8lf \t", rho_dust[l][Num_out_particle], vel_dust[l][Num_out_particle]);
			}
			fprintf(out_particle_file, "\n");
		}
	} while (Tm < T_end);

	fclose(out_particle_file);
	//	getchar();
	return 0;
}


