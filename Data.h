
// ����������������� ��������
double Rm, Zm, R0, c0, c_gas, Gam_liq, p_left, rho_left, p_right, rho_right, rho_dust_left, rho_dust_right, B, tau, Tm, Cg_0;
	double Xmax, Xmin;
	double rho_dust_0, mas_gas_0;
	double R_dust[10], mas_dustgas[10];
	double *p_s, *p_gas, *p_g; // ��������
	double *u_s, *u_p; // ��������
	double *rho_s, *rho_gas, *rho_p1; // ���������
	double h0, delta, rho_init;

	double *Cg; // ������������ 
	double *Beta, *S, *SS; // ������, �������� � ��������� ��������
	double Gam_g;
	long Rcount, Rn, RIm;
	const double p0=1.0e5;
	const double rho0_liq = 1000;
	const double rho0_gas = 1.14;
	const double v0 = 10;
	const double PI = 3.141592653589793238463;

	int Variant_type, Num_gas_particle_0, Type_of_state;
	
// ��������� ��������

	int Im1_gas, Im2_gas, Im1_dust, Im2_dust; // ����� ������
	int Num_dust_particle, Num_gas_particle, Num_dust_sort;
	int Jm; // ����� �����
	double *x_gas;// ���������� �������
	double *h_gas; // ������ �����������
	double *mas_gas; // ����� �������
	double *acc_gas; // ��������� �������
	double *vel_gas; // �������� �������
	double *e_gas;
	double *Rtau; // ��� �� ������� ��� ���������
	double T_end, T_out; // ����� ����� �� �������
	int Il, Ir; // ������ ��������� ������
	int File_int; // �������� ������ � ����
	double Geps; // �������� ���� �������
	double alpha, beta, eps, Cnu; // ��������� ����. ��������
	int *ind_gas; // ������ �������

	double *x_dust[10];// ���������� �������
	double *h_dust[10]; // ������ �����������
	double *rho_dust[10]; // ������ �����������
	double *mas_dust[10]; // ����� �������
	double *acc_dust[10]; // ��������� �������
	double *vel_dust[10]; // �������� �������
	double *t_stop[10];
	int *ind_dust[10]; // ������ �������

	double t_stop_init[10], a_rho[10], b_rho[10], a_vel[10], b_vel[10], a_vel_gas, b_vel_gas, a_rho_gas, b_rho_gas;


// ������� ��� GPU

	
    double *dev_P = 0;
	double *dev_rho = 0;
    double *dev_Pg = 0;	
    double *dev_R = 0;		
	double *R_new = 0;
    double *dev_S = 0;
    double *dev_SS = 0;
    double *dev_Cg = 0;	
    double *dev_Rtau = 0;	
	double *dev_RTm = 0;	

// ������

	double RSS;
	double R1, R2, R3, R4, R5, S1, S2, S3, S4, S5, SS1, SS2, SS3, SS4, SS5, Dl, Dl1, Dl2, Reps;


// Cell
	int Number_of_cell;
	double *x_cell;
	double *v_g_average, *v_d_average[10];
	double *rho_g_average, *rho_d_average[10], *eps_cell[10];
	double *t_stop_average[10], *e_g_average;
	double *x_av_new[10], *y_av_new, *x_av[10], *y_av, *Psi_av_new, *u_av_new[10], *v_av_new, *b_cell[10], *beta_cell;
	int *g_average_count, *d_average_count[10];
	double cell_width;

	int Flag_of_method;