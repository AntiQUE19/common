
    #include "stdio.h"
    #include "stdlib.h"
    #include "math.h"

//  The program minimises functional using maximum principle of Pontryagin.
//  Equalities are integrated by Runge-Kutta method version DP7-8
//  To solve linear equalities system Gauss method is used.

    //functions of ODE system and functional
	double f1(double, double, double, double, double);
	double f2(double, double, double, double, double);
	double f3(double, double, double, double, double);
	double f4(double, double, double, double, double);
	double B0(double, double, double, double, double);

	double max(double, double); //finding maximum of 2 numbers
	double min(double, double); //finding minimum of 2 numbers
	void gauss(double [2][3]); // Gauss method for linear equalities
	void runge(double, double, double, double); //RK method
	void runge_print(double, double, double, double);//RK method, result is written to file
	void back(double, double, double, double, double);// RK method: calculating backwards
	void runge_numbers(double *);// calculating Runge numbers in order to control solution
	void runge4(double, double, double, double);//RK for f1, f2, f3, f4

	double x, y, px, py, t; //variables on the present step
	double u;
    double xo, yo, pxo, pyo, bo;//initial values
    double end;//end of period of time
	double xT, yT, pxT, pyT;//auxilary variables for RKM
	double x_prev, y_prev, px_prev, py_prev, t_prev, bo_prev; //variables on the previous step
	double vec[] = {0.0, 0.01, 0.5, 1.5, 10.5}; //solution parameters
	double alpha;
	int num_alpha;

    //variables for RKM
	double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
	double q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13;
	double m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13;
	double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13;
	double l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13;

	double r1, r2, r3, r4, r5, h;
	double eps1, eps2;
    int n;
	int m;
	int q;
	int cntr;

	double pi;
	double precision;
    double ep;

	double a[2][3];
	double a11, a12, a21, a22;

	double px_n, py_n, x_n, y_n, r_bo;

	int i;

	double ge0, ge1;
	double t_ch;

	double err, eps, fac, facmin, facmax, alpha;//parameters of solution
    double Rx, Ry, Rpx, Rpy; //Runge numbers
    int j;//variable for solution
    double L, ls, d, tmp;//extra variables

    //coefficients for RKM
	double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13;/* alpha i */
    double b21, b31, b41, b51, b61, b71, b81, b91, b101, b111, b121, b131;/*beta i(2to13) j=1 */
    double b32, b42, b52, b62, b72, b82, b92, b102, b112, b122, b132; /*beta i (3to13) j=2 */
    double b43, b53, b63, b73, b83, b93, b103, b113, b123, b133; /* beta i (4to13) j=3 */
    double b54, b64, b74, b84, b94, b104, b114, b124, b134; /* beta i (5to13) j=4 */
    double b65, b75, b85, b95, b105, b115, b125, b135; /* beta i (6to13) j=5 */
    double b76, b86, b96, b106, b116, b126, b136; /*beta i (7to13) j=6 */
    double b87, b97, b107, b117, b127, b137; /* beta i (8 to13) j=7 */
    double b98, b108, b118, b128, b138; /*beta i (9to13) j=8 */
    double b109, b119, b129, b139; /* beta i (10to13) j=9 */
    double b1110, b1210, b1310; /* beta i (11to13) j=10 */
    double b1211, b1311, b1312;/* beta i (12,13) j=11 */
    double g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13; /* gammas */
    double gt1, gt2, gt3, gt4, gt5, gt6, gt7, gt8, gt9, gt10, gt11, gt12, gt13; /* gammas ~ */

int main(void)
{
	FILE *out;
	double *b;
    //memory allocation for variables
	if (!(b = (double*)malloc(12*sizeof(double))))
    {
           printf("Memory allocation failure\n");
           return (-1);
    }
    //initiation of cooefficients of RKM
    a1=0.0, a2=1.0/18.0, a3=1.0/12.0, a4=1.0/8.0, a5=5.0/16.0, a6=3.0/8.0, a7=59.0/400.0, a8=93.0/200.0;
    a9=5490023248.0/9719169821.0, a10=13.0/20.0, a11=1201146811.0/1299019798.0, a12=1.0, a13=1.0;/* alpha i */
    b21=1.0/18.0, b31=1.0/48.0, b41=1.0/32.0, b51=5.0/16.0, b61=3.0/80.0, b71=29443841.0/614563906.0;
    b81=16016141.0/946692911.0, b91=39632708.0/573591083.0, b101=246121993.0/1340847787.0;
    b111=-1028468189.0/846180014.0, b121=185892177.0/718116043.0, b131=403863854.0/491063109.0;/*beta i(2to13) j=1 */
    b32=1.0/16.0, b42=0.0, b52=0.0, b62=0.0, b72=0.0, b82=0.0, b92=0.0, b102=0.0, b112=0.0;
    b122=0.0, b132=0.0; /*beta i (3to13) j=2 */
    b43=3.0/32.0, b53=-75.0/64.0, b63=0.0, b73=0.0, b83=0.0, b93=0.0, b103=0.0, b113=0.0;
    b123=0.0, b133=0.0; /* beta i (4to13) j=3 */
    b54=75.0/64.0, b64=3.0/16.0, b74=77736538.0/692538347.0, b84=61564180.0/158732637.0, b94=-433636366.0/683701615.0;
    b104=-37695042795.0/15268766246.0, b114=8478235783.0/508512852.0, b124=-3185094517.0/667107341.0;
    b134=-5068492393.0/434740067; /* beta i (5to13) j=4 */
    b65=3.0/20.0, b75=-28693883.0/1125000000.0, b85=22789713.0/633445777.0, b95=-421739975.0/2616292301.0;
    b105=-309121744.0/1061227803.0, b115=1311729495.0/1432422823.0, b125=-477755414.0/1098053517.0;
    b135=-411421997.0/543043805.0; /* beta i (6to13) j=5 */
    b76=23124283.0/1800000000.0, b86=545815736.0/2771057229.0, b96=100302831.0/723423059.0;
    b106=-12992083.0/490766935.0, b116=-10304129995.0/1701304382.0, b126=-703635378.0/230739211.0;
    b136=652783627.0/914296604.0; /*beta i (7to13) j=6 */
    b87=-180193667.0/1043307555.0, b97=790204164.0/839813087.0, b107=6005943493.0/2108947869.0;
    b117=-48777925059.0/3047939560.0, b127=5731566787.0/1027545527.0, b137=11173962825.0/925320556.0; /* beta i (8 to13) j=7 */
    b98=800635310.0/3783071287.0, b108=393006217.0/1396673457.0, b118=15336726248.0/1032824649.0;
    b128=5232866602.0/850066563.0, b138=-13158990841.0/6184727034.0; /*beta i (9to13) j=8 */
    b109=123872331.0/1001029789.0, b119=-45442868181.0/3398467696.0, b129=-4093664535.0/808688257.0;
    b139=3936647629.0/1978049680.0; /* beta i (10to13) j=9 */
    b1110=3065993473.0/597172653.0, b1210=3962137247.0/1805957418.0, b1310=-160528059.0/685178525.0; /* beta i (11to13) j=10 */
    b1211=65686358.0/487910083.0, b1311=248638103.0/1413531060.0, b1312=0.0;/* beta i (12,13) j=11 */
    g1=14005451.0/335480064.0, g2=0.0, g3=0.0, g4=0.0, g5=0.0, g6=-59238493.0/1068277825.0;
    g7=181606767.0/758867731.0, g8=561292985.0/797845732.0, g9=-1041891430.0/1371343529.0, g10=760417239.0/1151165299.0;
    g11=118820643.0/751138087.0, g12=-528747749.0/2220607170.0, g13=1.0/4.0; /* gammas */
    gt1=13451932.0/455176623.0, gt2=0.0, gt3=0.0, gt4=0.0, gt5=0.0, gt6=-808719846.0/976000145.0;
    gt7=1757004468.0/5645159321.0, gt8=656045339.0/265891186.0, gt9=-3867574721.0/1518517206.0;
    gt10=465885868.0/322736535.0, gt11=53011238.0/667516719.0, gt12=2.0/45.0, gt13=0.0; /* gammas ~ */

	xo = 0.0;
    t = 0.0;
	alpha = 0;
	end = M_PI/2;
	facmin = 0.5;
    facmax = 2.5;
    fac = 0.85;
    eps = 1.e-12;
    precision = 1.e-8;
    ep = 1.e-14;

	out = fopen("out.txt" ,"w");
	for (i = 0; i < 5; i++)
    {
		alpha = vec[i];
        if (alpha == 0)
        {
			pxo = 0.0;
			pyo = 4/ (M_PI - 4);
		}
        yo = -pyo;
        n = 0;
		m = 0;
		q = 0;

        runge(xo, yo, pxo, pyo);

		px_n = px;
		py_n = py;
		x_n = x;
		y_n = y;

		while ((fabs(x-1.0) > precision) || (fabs(py) > precision))
        {
			runge(xo, yo, pxo+ep, pyo);
			a[0][0] = (px - px_n)/ep;
			a[1][0] = (py - py_n)/ep;

			runge(xo, yo, pxo, pyo+ep);
			a[0][1] = (px - px_n)/ep;
			a[1][1] = (py - py_n)/ep;

			a[0][2] = -(x_n-1.0) + a[0][0] * pxo + a[0][1] * pyo;
			a[1][2] = -(py_n) + a[1][0] * pxo + a[1][1] * pyo;

			gauss(a);

			pxo = a[0][2];
            pyo = a[1][2];

            runge(xo, yo, pxo, pyo);

			px_n = px;
			py_n = py;
			x_n = x;
			y_n = y;
            m++;
		}

        fprintf(out, "\n\n");
		fprintf(out, "\n\t * * * alpha = %lf * * *\n\n", alpha);
		runge_print(xo, yo, pxo, pyo);
		fprintf(out, "\n\tStarting and shooting parameters:\n\n x0 = %.15lf,\n y0 = %.15lf,\n px0 = %.15lf,\n py0 = %.15lf\n", xo, yo, pxo, pyo);
		fprintf(out, "\n\tIt took %d steps of shooting \n", m);
		fprintf(out, "\n t = %.20lf\n x0 = %.20lf\n y0 = %.20lf\n px0 = %.20lf\n py0 = %.15lf\n",0., xo, yo, pxo, pyo);
		fprintf(out, "\n t = %.20lf\n x = %.20lf\n y = %.20lf\n px = %.20lf\n py = %.15lf\n Bo = %.15lf\n", t, x, y, px, py, bo);
		fprintf(out, "\n\tPrecision vector: \n\n |x - 1.0| = %5.20lf\n |py| = %5.20lf\n\n", fabs(x - 1.0), fabs(py));
		fprintf(out, " \n\t Calculating Backwards: \n");
		back(x, y, px, py, t);
		fprintf(out, "\n\tAnswer is: \n\n t = %5.20lf\n x = %5.20lf\n y = %5.20lf\n px = %5.20lf\n py = %5.20lf\n\n", t, x, y, px, py);
		fprintf(out, "\n\t Should be: \n\n t = %5.20lf\n x = %5.20lf\n y = %5.20lf\n px = %5.20lf\n py = %5.20lf\n\n", 0.0, xo, yo, pxo, pyo);
		fprintf(out, "\n\tPrecision vector:\n\n|x - xo| = %5.20lf\n |y - yo| =  %5.20lf\n |px - pxo| = %5.20lf\n |py - pyo| =%5.20lf", x - xo, y - yo, px - pxo, py - pyo);
    }

	fclose(out);
	runge_numbers(b);
	printf("\nSee result in 'out.txt'\n\n");
	free(b);
	return (0);
}

double max(double a, double b)
{
	return a > b ? a : b;
}

double min(double a, double b)
{
	return a < b ? a : b;
}

double f1(double t,double x,double y, double px, double py)
{
	double f;
	t = t;
	x = x;
	px = px;
	py = py;

	f = y;
	return (f);
}

double f2(double t,double x,double y, double px, double py)
{

	double f;
	t = t;
	x = x;
	y = y;
	px = px;
	py = py;

    f = py-x/(1+alpha*t*t);
	return (f);
}

double f3(double t,double x,double y, double px, double py)
{
	double f;
	t = t;
	x = x;
	y = y;
	px = px;
	py = py;

    f = py/(1+alpha*t*t);
	return (f);
}

double f4(double t, double x, double y, double px, double py)
{
	double f;
	t = t;
	x = x;
	y = y;
	px = px;
	py = py;

	f = - px;
	return (f);
}

double B0(double t,double x,double y, double px, double py)
{
	double f;
	t = t;
	x = x;
	y = y;
	px = px;
	py = py;

	f = py * py;
	return (f);
}

void runge(double r_xo, double r_yo, double r_pxo, double r_pyo)
{
	x_prev = r_xo;
	y_prev = r_yo;
	px_prev = r_pxo;
	py_prev = r_pyo;

	t_prev = 0.0;
	h = 0.01;
	t = 0.0;
	q = 0;
    cntr = 1;

    x = r_xo;
    y = r_yo;
    px = r_pxo;
    py = r_pyo;

    while(t < end)
    {
        n++;
        k1=f1(t+a1*h, x, y, px, py);
        q1=f2(t+a1*h, x, y, px, py);
        m1=f3(t+a1*h, x, y, px, py);
        p1=f4(t+a1*h, x, y, px, py);

        k2=f1(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        q2=f2(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        m2=f3(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        p2=f4(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);

        k3=f1(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        q3=f2(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        m3=f3(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        p3=f4(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));

        k4=f1(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        q4=f2(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        m4=f3(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        p4=f4(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));

        k5=f1(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        q5=f2(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        m5=f3(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        p5=f4(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));

        k6=f1(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        q6=f2(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        m6=f3(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        p6=f4(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));

        k7=f1(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        q7=f2(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        m7=f3(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        p7=f4(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));

        k8=f1(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        q8=f2(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        m8=f3(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        p8=f4(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));

        k9=f1(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        q9=f2(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        m9=f3(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        p9=f4(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));

        k10=f1(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        q10=f2(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        m10=f3(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        p10=f4(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));

        k11=f1(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        q11=f2(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        m11=f3(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        p11=f4(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));

        k12=f1(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        q12=f2(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        m12=f3(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        p12=f4(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));

        k13=f1(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        q13=f2(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        m13=f3(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        p13=f4(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));

		r1 = fabs(h*(k1*(g1-gt1)+k2*(g2-gt2)+k3*(g3-gt3)+k4*(g4-gt4)+k5*(g5-gt5)+k6*(g6-gt6)+k7*(g7-gt7)+k8*(g8-gt8)+k9*(g9-gt9)+k10*(g10-gt10)+k11*(g11-gt11)+k12*(g12-gt12)+k13*(g13-gt13)));
		r2 = fabs(h*(q1*(g1-gt1)+q2*(g2-gt2)+q3*(g3-gt3)+q4*(g4-gt4)+q5*(g5-gt5)+q6*(g6-gt6)+q7*(g7-gt7)+q8*(g8-gt8)+q9*(g9-gt9)+q10*(g10-gt10)+q11*(g11-gt11)+q12*(g12-gt12)+q13*(g13-gt13)));
		r3 = fabs(h*(m1*(g1-gt1)+m2*(g2-gt2)+m3*(g3-gt3)+m4*(g4-gt4)+m5*(g5-gt5)+m6*(g6-gt6)+m7*(g7-gt7)+m8*(g8-gt8)+m9*(g9-gt9)+m10*(g10-gt10)+m11*(g11-gt11)+m12*(g12-gt12)+m13*(g13-gt13)));
		r4 = fabs(h*(p1*(g1-gt1)+p2*(g2-gt2)+p3*(g3-gt3)+p4*(g4-gt4)+p5*(g5-gt5)+p6*(g6-gt6)+p7*(g7-gt7)+p8*(g8-gt8)+p9*(g9-gt9)+p10*(g10-gt10)+p11*(g11-gt11)+p12*(g12-gt12)+p13*(g13-gt13)));

		x += h*(g1*k1+g2*k2+g3*k3+g4*k4+g5*k5+g6*k6+g7*k7+g8*k8+g9*k9+g10*k10+g11*k11+g12*k12+g13*k13);;
		y += h*(g1*q1+g2*q2+g3*q3+g4*q4+g5*q5+g6*q6+g7*q7+g8*q8+g9*q9+g10*q10+g11*q11+g12*q12+g13*q13);;
		px += h*(g1*m1+g2*m2+g3*m3+g4*m4+g5*m5+g6*m6+g7*m7+g8*m8+g9*m9+g10*m10+g11*m11+g12*m12+g13*m13);;
		py += h*(g1*p1+g2*p2+g3*p3+g4*p4+g5*p5+g6*p6+g7*p7+g8*p8+g9*p9+g10*p10+g11*p11+g12*p12+g13*p13);;

        err = max( max(r1, r2), max(r3,r4));
        if(err > eps)
        {
            h = h * min( max( facmin, fac*pow(eps/err, 1/9.0)), facmax);

            t = t_prev;
            x = x_prev;
            y = y_prev;
            px = px_prev;
			py = py_prev;
            n--;
        }
        else
        {
            x_prev = x;
            y_prev = y;
            px_prev = px;
            py_prev = py;
            t_prev = t;
            t+=h;
        }

		if(t + h == end)
        {
            cntr=0;
        }
        if((t + h > end)&&(cntr==1))
        {
            h = end - t;
            cntr = 0;
        }
    }
	return;
}

void runge_print(double r_xo, double r_yo, double r_pxo, double r_pyo)
{

    FILE *resultx, *resulty, *resultpx, *resultpy, *resultbo;

	resultx = fopen("resultx.txt" ,"w");
	resulty = fopen("resulty.txt" ,"w");
	resultpx = fopen("resultpx.txt" ,"w");
	resultpy = fopen("resultpy.txt" ,"w");
	resultbo = fopen("resultbo.txt", "w");

	x_prev = r_xo;
	y_prev = r_yo;
	px_prev = r_pxo;
	py_prev = r_pyo;
	bo_prev = r_bo;

	t_prev = 0;
	h = 0.01;
	t = 0;
	q = 0;
	end = M_PI/2;
    cntr = 1;
    tmp = 0;

    x = r_xo;
    y = r_yo;
    px = r_pxo;
    py = r_pyo;
    bo = bo_prev;

    while(t < end)
    {

        n++;
        k1=f1(t+a1*h, x, y, px, py);
        q1=f2(t+a1*h, x, y, px, py);
        m1=f3(t+a1*h, x, y, px, py);
        p1=f4(t+a1*h, x, y, px, py);
        l1=B0(t+a1*h, x, y, px, py);

        k2=f1(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        q2=f2(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        m2=f3(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        p2=f4(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        l2=B0(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);

        k3=f1(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        q3=f2(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        m3=f3(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        p3=f4(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        l3=B0(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));

        k4=f1(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        q4=f2(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        m4=f3(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        p4=f4(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        l4=B0(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));

        k5=f1(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        q5=f2(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        m5=f3(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        p5=f4(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        l5=B0(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));

        k6=f1(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        q6=f2(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        m6=f3(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        p6=f4(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        l6=B0(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));

        k7=f1(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        q7=f2(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        m7=f3(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        p7=f4(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        l7=B0(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));

        k8=f1(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        q8=f2(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        m8=f3(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        p8=f4(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        l8=B0(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));

        k9=f1(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        q9=f2(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        m9=f3(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        p9=f4(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        l9=B0(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));

        k10=f1(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        q10=f2(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        m10=f3(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        p10=f4(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        l10=B0(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));

        k11=f1(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        q11=f2(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        m11=f3(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        p11=f4(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        l11=B0(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));

        k12=f1(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        q12=f2(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        m12=f3(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        p12=f4(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        l12=B0(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));

        k13=f1(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        q13=f2(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        m13=f3(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        p13=f4(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        l13=B0(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));

		r1 = fabs(h*(k1*(g1-gt1)+k2*(g2-gt2)+k3*(g3-gt3)+k4*(g4-gt4)+k5*(g5-gt5)+k6*(g6-gt6)+k7*(g7-gt7)+k8*(g8-gt8)+k9*(g9-gt9)+k10*(g10-gt10)+k11*(g11-gt11)+k12*(g12-gt12)+k13*(g13-gt13)));
		r2 = fabs(h*(q1*(g1-gt1)+q2*(g2-gt2)+q3*(g3-gt3)+q4*(g4-gt4)+q5*(g5-gt5)+q6*(g6-gt6)+q7*(g7-gt7)+q8*(g8-gt8)+q9*(g9-gt9)+q10*(g10-gt10)+q11*(g11-gt11)+q12*(g12-gt12)+q13*(g13-gt13)));
		r3 = fabs(h*(m1*(g1-gt1)+m2*(g2-gt2)+m3*(g3-gt3)+m4*(g4-gt4)+m5*(g5-gt5)+m6*(g6-gt6)+m7*(g7-gt7)+m8*(g8-gt8)+m9*(g9-gt9)+m10*(g10-gt10)+m11*(g11-gt11)+m12*(g12-gt12)+m13*(g13-gt13)));
		r4 = fabs(h*(p1*(g1-gt1)+p2*(g2-gt2)+p3*(g3-gt3)+p4*(g4-gt4)+p5*(g5-gt5)+p6*(g6-gt6)+p7*(g7-gt7)+p8*(g8-gt8)+p9*(g9-gt9)+p10*(g10-gt10)+p11*(g11-gt11)+p12*(g12-gt12)+p13*(g13-gt13)));
        r5 = fabs(h*(l1*(g1-gt1)+l2*(g2-gt2)+l3*(g3-gt3)+l4*(g4-gt4)+l5*(g5-gt5)+l6*(g6-gt6)+l7*(g7-gt7)+l8*(g8-gt8)+l9*(g9-gt9)+l10*(g10-gt10)+l11*(g11-gt11)+l12*(g12-gt12)+l13*(g13-gt13)));

		x += h*(g1*k1+g2*k2+g3*k3+g4*k4+g5*k5+g6*k6+g7*k7+g8*k8+g9*k9+g10*k10+g11*k11+g12*k12+g13*k13);;
		y += h*(g1*q1+g2*q2+g3*q3+g4*q4+g5*q5+g6*q6+g7*q7+g8*q8+g9*q9+g10*q10+g11*q11+g12*q12+g13*q13);;
		px += h*(g1*m1+g2*m2+g3*m3+g4*m4+g5*m5+g6*m6+g7*m7+g8*m8+g9*m9+g10*m10+g11*m11+g12*m12+g13*m13);;
		py += h*(g1*p1+g2*p2+g3*p3+g4*p4+g5*p5+g6*p6+g7*p7+g8*p8+g9*p9+g10*p10+g11*p11+g12*p12+g13*p13);;
        bo += h*(g1*l1+g2*l2+g3*l3+g4*l4+g5*l5+g6*l6+g7*l7+g8*l8+g9*l9+g10*l10+g11*l11+g12*l12+g13*l13);;

        err = max( max(r1, r2), max(r3,r4));

        if(err>eps)
        {
            h = h * min( max( facmin, fac*pow(eps/err, 1/9.0)), facmax);

            t = t_prev;
            x = x_prev;
            y = y_prev;
            px = px_prev;
			py = py_prev;
			bo = bo_prev;
            n--;
        }
        else
        {
            x_prev = x;
            y_prev = y;
            px_prev = px;
            py_prev = py;
            bo_prev = bo;
            t_prev = t;
            t+=h;

        }

		if(t + h == end)
        {
            cntr = 0;
        }
        if((t + h > end)&&(cntr == 1))
        {
            h = end - t;
            cntr = 0;
        }

		fprintf(resultx,"%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", t, x, y, px, py, bo);
		fprintf(resulty,"%lf %lf\n", t, y);
		fprintf(resultpx,"%lf %lf\n", t, px);
		fprintf(resultpy,"%lf %lf\n", t, py);
		fprintf(resultbo,"%lf %lf\n", t, bo);

	}
	return;
}

void gauss(double a[2][3])
{
	int j, ii, jj;
	double buf[3];
	double epsil = 1.e-90;

	if (fabs(a[0][0]) < epsil)
	{
		if (fabs(a[1][0]) < epsil)
        {
			for (ii = 0; ii < 2; ii++)
            {
				//printf("\t");
					for (jj = 0; jj < 3; jj++)
                    {
						//printf ("%e ",a[ii][jj]);
					}
				//printf("\n");
			}
			return ;
		}
		else
        {
			for (jj = 0; jj < 3; jj++)
			{
				buf[jj] = a[1][jj];
				a[1][jj] = a[0][jj];
				a[0][jj] = buf[jj];
			}
		}
	}

	for (j = 2; j >= 0; j--)
    {
		a[0][j] /= a[0][0];
	}

	for (j = 2; j >= 0 ; j--)
    {
		a[1][j] -= a[0][j] * a[1][0];
	}

	if (fabs(a[1][1])<epsil)
    {
		for (ii = 0; ii < 2; ii++)
        {
			//printf("\t");
				for (jj = 0; jj < 3; jj++)
				{
					//printf ("%e ",a[ii][jj]);
				}
			//printf("\n");
		}
		return;
	}

	if (fabs(a[1][1]) > epsil)
    {
		for (j = 2; j > 0; j--)
		{
			a[1][j] /= a[1][1];
		}
	}

	for (j = 2; j >= 1; j--)
    {
		a[0][j] -= a[1][j] * a[0][1];
	}

	return;
}

void back(double r_xo, double r_yo, double r_pxo, double r_pyo, double r_end)
{
	x_prev = r_xo;
	y_prev = r_yo;
	px_prev = r_pxo;
	py_prev = r_pyo;
	t_prev = r_end;

	h = -0.0001;

	cntr = 1;

    x = r_xo;
    y = r_yo;
    px = r_pxo;
    py = r_pyo;
    t = r_end;

    while(t > 0)
    {
        n++;
        k1=f1(t+a1*h, x, y, px, py);
        q1=f2(t+a1*h, x, y, px, py);
        m1=f3(t+a1*h, x, y, px, py);
        p1=f4(t+a1*h, x, y, px, py);
        //l1=B0(t+a1*h, x, y, px, py);

        k2=f1(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        q2=f2(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        m2=f3(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        p2=f4(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
       // l2=B0(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);

        k3=f1(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        q3=f2(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        m3=f3(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        p3=f4(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        //l3=B0(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));

        k4=f1(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        q4=f2(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        m4=f3(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        p4=f4(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
       // l4=B0(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));

        k5=f1(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        q5=f2(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        m5=f3(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        p5=f4(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        //l5=B0(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));

        k6=f1(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        q6=f2(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        m6=f3(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        p6=f4(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        //l6=B0(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));

        k7=f1(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        q7=f2(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        m7=f3(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        p7=f4(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
       // l7=B0(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));

        k8=f1(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        q8=f2(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        m8=f3(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        p8=f4(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
       //l8=B0(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));

        k9=f1(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        q9=f2(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        m9=f3(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        p9=f4(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));

        k10=f1(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        q10=f2(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        m10=f3(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        p10=f4(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
       // l10=B0(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));

        k11=f1(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        q11=f2(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        m11=f3(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        p11=f4(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        //l11=B0(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));

        k12=f1(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        q12=f2(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        m12=f3(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        p12=f4(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
       // l12=B0(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));

        k13=f1(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        q13=f2(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        m13=f3(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        p13=f4(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
       // l13=B0(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));

		r1 = fabs(h*(k1*(g1-gt1)+k2*(g2-gt2)+k3*(g3-gt3)+k4*(g4-gt4)+k5*(g5-gt5)+k6*(g6-gt6)+k7*(g7-gt7)+k8*(g8-gt8)+k9*(g9-gt9)+k10*(g10-gt10)+k11*(g11-gt11)+k12*(g12-gt12)+k13*(g13-gt13)));
		r2 = fabs(h*(q1*(g1-gt1)+q2*(g2-gt2)+q3*(g3-gt3)+q4*(g4-gt4)+q5*(g5-gt5)+q6*(g6-gt6)+q7*(g7-gt7)+q8*(g8-gt8)+q9*(g9-gt9)+q10*(g10-gt10)+q11*(g11-gt11)+q12*(g12-gt12)+q13*(g13-gt13)));
		r3 = fabs(h*(m1*(g1-gt1)+m2*(g2-gt2)+m3*(g3-gt3)+m4*(g4-gt4)+m5*(g5-gt5)+m6*(g6-gt6)+m7*(g7-gt7)+m8*(g8-gt8)+m9*(g9-gt9)+m10*(g10-gt10)+m11*(g11-gt11)+m12*(g12-gt12)+m13*(g13-gt13)));
		r4 = fabs(h*(p1*(g1-gt1)+p2*(g2-gt2)+p3*(g3-gt3)+p4*(g4-gt4)+p5*(g5-gt5)+p6*(g6-gt6)+p7*(g7-gt7)+p8*(g8-gt8)+p9*(g9-gt9)+p10*(g10-gt10)+p11*(g11-gt11)+p12*(g12-gt12)+p13*(g13-gt13)));
       // r5 = fabs(h*(l1*(g1-gt1)+l2*(g2-gt2)+l3*(g3-gt3)+l4*(g4-gt4)+l5*(g5-gt5)+l6*(g6-gt6)+l7*(g7-gt7)+l8*(g8-gt8)+l9*(g9-gt9)+l10*(g10-gt10)+l11*(g11-gt11)+l12*(g12-gt12)+l13*(g13-gt13)));

		x += h*(g1*k1+g2*k2+g3*k3+g4*k4+g5*k5+g6*k6+g7*k7+g8*k8+g9*k9+g10*k10+g11*k11+g12*k12+g13*k13);;
		y += h*(g1*q1+g2*q2+g3*q3+g4*q4+g5*q5+g6*q6+g7*q7+g8*q8+g9*q9+g10*q10+g11*q11+g12*q12+g13*q13);;
		px += h*(g1*m1+g2*m2+g3*m3+g4*m4+g5*m5+g6*m6+g7*m7+g8*m8+g9*m9+g10*m10+g11*m11+g12*m12+g13*m13);;
		py += h*(g1*p1+g2*p2+g3*p3+g4*p4+g5*p5+g6*p6+g7*p7+g8*p8+g9*p9+g10*p10+g11*p11+g12*p12+g13*p13);;
        //bo += h*(g1*l1+g2*l2+g3*l3+g4*l4+g5*l5+g6*l6+g7*l7+g8*l8+g9*l9+g10*l10+g11*l11+g12*l12+g13*l13);;

        err = max( max(r1, r2), max(r3,r4));

        if(err > eps)
        {
            h = h * min( max( facmin, fac*pow(eps/err, 1/9.0)), facmax);
            t = t_prev;
            x = x_prev;
            y = y_prev;
            px = px_prev;
			py = py_prev;
            n--;
        }
        else
        {
            x_prev = x;
            y_prev = y;
            px_prev = px;
            py_prev = py;
            t_prev = t;
            t+=h;
        }
		if(t + h == 0.0)
        {
            cntr=0;
        }
        if((t + h < 0)&&(cntr==1))
        {
            h = - t;
            cntr = 0;
        }
    }
	return;
}

void runge4(double r_xo, double r_yo, double r_pxo, double r_pyo)
{
	x_prev = r_xo;
	y_prev = r_yo;
	px_prev = r_pxo;
	py_prev = r_pyo;

	t_prev = 0.0;
	h = 0.01;
	t = 0.0;
	q = 0;
    cntr = 1;

    x = r_xo;
    y = r_yo;
    px = r_pxo;
    py = r_pyo;

    while(t < end)
    {
        n++;
        k1=f1(t+a1*h, x, y, px, py);
        q1=f2(t+a1*h, x, y, px, py);
        m1=f3(t+a1*h, x, y, px, py);
        p1=f4(t+a1*h, x, y, px, py);

        k2=f1(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        q2=f2(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        m2=f3(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);
        p2=f4(t+a2*h, x+h*b21*k1, y+h*b21*q1, px+h*b21*m1, py+h*b21*p1);

        k3=f1(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        q3=f2(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        m3=f3(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));
        p3=f4(t+a3*h, x+h*(b31*k1+b32*k2), y+h*(b31*q1+b32*q2), px+h*(b31*m1+b32*m2), py+h*(b31*p1+b32*p2));

        k4=f1(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        q4=f2(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        m4=f3(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));
        p4=f4(t+a4*h, x+h*(b41*k1+b42*k2+b43*k3), y+h*(b41*q1+b42*q2+b43*q3), px+h*(b41*m1+b42*m2+b43*m3), py+h*(b41*p1+b42*p2+b43*p3));

        k5=f1(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        q5=f2(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        m5=f3(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));
        p5=f4(t+a5*h, x+h*(b51*k1+b52*k2+b53*k3+b54*k4), y+h*(b51*q1+b52*q2+b53*q3+b54*q4), px+h*(b51*m1+b52*m2+b53*m3+b54*m4), py+h*(b51*p1+b52*p2+b53*p3+b54*p4));

        k6=f1(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        q6=f2(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        m6=f3(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));
        p6=f4(t+a6*h, x+h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5), y+h*(b61*q1+b62*q2+b63*q3+b64*q4+b65*q5), px+h*(b61*m1+b62*m2+b63*m3+b64*m4+b65*m5), py+h*(b61*p1+b62*p2+b63*p3+b64*p4+b65*p5));

        k7=f1(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        q7=f2(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        m7=f3(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));
        p7=f4(t+a7*h, x+h*(b71*k1+b72*k2+b73*k3+b74*k4+b75*k5+b76*k6), y+h*(b71*q1+b72*q2+b73*q3+b74*q4+b75*q5+b76*q6), px+h*(b71*m1+b72*m2+b73*m3+b74*m4+b75*m5+b76*m6), py+h*(b71*p1+b72*p2+b73*p3+b74*p4+b75*p5+b76*p6));

        k8=f1(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        q8=f2(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        m8=f3(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));
        p8=f4(t+a8*h, x+h*(b81*k1+b82*k2+b83*k3+b84*k4+b85*k5+b86*k6+b87*k7), y+h*(b81*q1+b82*q2+b83*q3+b84*q4+b85*q5+b86*q6+b87*q7), px+h*(b81*m1+b82*m2+b83*m3+b84*m4+b85*m5+b86*m6+b87*m7), py+h*(b81*p1+b82*p2+b83*p3+b84*p4+b85*p5+b86*p6+b87*p7));

        k9=f1(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        q9=f2(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        m9=f3(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));
        p9=f4(t+a9*h, x+h*(b91*k1+b92*k2+b93*k3+b94*k4+b95*k5+b96*k6+b97*k7+b98*k8), y+h*(b91*q1+b92*q2+b93*q3+b94*q4+b95*q5+b96*q6+b97*q7+b98*q8), px+h*(b91*m1+b92*m2+b93*m3+b94*m4+b95*m5+b96*m6+b97*m7+b98*m8), py+h*(b91*p1+b92*p2+b93*p3+b94*p4+b95*p5+b96*p6+b97*p7+b98*p8));

        k10=f1(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        q10=f2(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        m10=f3(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));
        p10=f4(t+a10*h, x+h*(b101*k1+b102*k2+b103*k3+b104*k4+b105*k5+b106*k6+b107*k7+b108*k8+b109*k9), y+h*(b101*q1+b102*q2+b103*q3+b104*q4+b105*q5+b106*q6+b107*q7+b108*q8+b109*q9), px+h*(b101*m1+b102*m2+b103*m3+b104*m4+b105*m5+b106*m6+b107*m7+b108*m8+b109*m9), py+h*(b101*p1+b102*p2+b103*p3+b104*p4+b105*p5+b106*p6+b107*p7+b108*p8+b109*p9));

        k11=f1(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        q11=f2(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        m11=f3(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));
        p11=f4(t+a11*h, x+h*(b111*k1+b112*k2+b113*k3+b114*k4+b115*k5+b116*k6+b117*k7+b118*k8+b119*k9+b1110*k10), y+h*(b111*q1+b112*q2+b113*q3+b114*q4+b115*q5+b116*q6+b117*q7+b118*q8+b119*q9+b1110*q10), px+h*(b111*m1+b112*m2+b113*m3+b114*m4+b115*m5+b116*m6+b117*m7+b118*m8+b119*m9+b1110*m10), py+h*(b111*p1+b112*p2+b113*p3+b114*p4+b115*p5+b116*p6+b117*p7+b118*p8+b119*p9+b1110*p10));

        k12=f1(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        q12=f2(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        m12=f3(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));
        p12=f4(t+a12*h, x+h*(b121*k1+b122*k2+b123*k3+b124*k4+b125*k5+b126*k6+b127*k7+b128*k8+b129*k9+b1210*k10+b1211*k11), y+h*(b121*q1+b122*q2+b123*q3+b124*q4+b125*q5+b126*q6+b127*q7+b128*q8+b129*q9+b1210*q10+b1211*q11), px+h*(b121*m1+b122*m2+b123*m3+b124*m4+b125*m5+b126*m6+b127*m7+b128*m8+b129*m9+b1210*m10+b1211*m11), py+h*(b121*p1+b122*p2+b123*p3+b124*p4+b125*p5+b126*p6+b127*p7+b128*p8+b129*p9+b1210*p10+b1211*p11));

        k13=f1(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        q13=f2(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        m13=f3(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));
        p13=f4(t+a13*h, x+h*(b131*k1+b132*k2+b133*k3+b134*k4+b135*k5+b136*k6+b137*k7+b138*k8+b139*k9+b1310*k10+b1311*k11+b1312*k12), y+h*(b131*q1+b132*q2+b133*q3+b134*q4+b135*q5+b136*q6+b137*q7+b138*q8+b139*q9+b1310*q10+b1311*q11+b1312*q12), px+h*(b131*m1+b132*m2+b133*m3+b134*m4+b135*m5+b136*m6+b137*m7+b138*m8+b139*m9+b1310*m10+b1311*m11+b1312*m12), py+h*(b131*p1+b132*p2+b133*p3+b134*p4+b135*p5+b136*p6+b137*p7+b138*p8+b139*p9+b1310*p10+b1311*p11+b1312*p12));

		r1 = fabs(h*(k1*(g1-gt1)+k2*(g2-gt2)+k3*(g3-gt3)+k4*(g4-gt4)+k5*(g5-gt5)+k6*(g6-gt6)+k7*(g7-gt7)+k8*(g8-gt8)+k9*(g9-gt9)+k10*(g10-gt10)+k11*(g11-gt11)+k12*(g12-gt12)+k13*(g13-gt13)));
		r2 = fabs(h*(q1*(g1-gt1)+q2*(g2-gt2)+q3*(g3-gt3)+q4*(g4-gt4)+q5*(g5-gt5)+q6*(g6-gt6)+q7*(g7-gt7)+q8*(g8-gt8)+q9*(g9-gt9)+q10*(g10-gt10)+q11*(g11-gt11)+q12*(g12-gt12)+q13*(g13-gt13)));
		r3 = fabs(h*(m1*(g1-gt1)+m2*(g2-gt2)+m3*(g3-gt3)+m4*(g4-gt4)+m5*(g5-gt5)+m6*(g6-gt6)+m7*(g7-gt7)+m8*(g8-gt8)+m9*(g9-gt9)+m10*(g10-gt10)+m11*(g11-gt11)+m12*(g12-gt12)+m13*(g13-gt13)));
		r4 = fabs(h*(p1*(g1-gt1)+p2*(g2-gt2)+p3*(g3-gt3)+p4*(g4-gt4)+p5*(g5-gt5)+p6*(g6-gt6)+p7*(g7-gt7)+p8*(g8-gt8)+p9*(g9-gt9)+p10*(g10-gt10)+p11*(g11-gt11)+p12*(g12-gt12)+p13*(g13-gt13)));

		x += h*(g1*k1+g2*k2+g3*k3+g4*k4+g5*k5+g6*k6+g7*k7+g8*k8+g9*k9+g10*k10+g11*k11+g12*k12+g13*k13);;
		y += h*(g1*q1+g2*q2+g3*q3+g4*q4+g5*q5+g6*q6+g7*q7+g8*q8+g9*q9+g10*q10+g11*q11+g12*q12+g13*q13);;
		px += h*(g1*m1+g2*m2+g3*m3+g4*m4+g5*m5+g6*m6+g7*m7+g8*m8+g9*m9+g10*m10+g11*m11+g12*m12+g13*m13);;
		py += h*(g1*p1+g2*p2+g3*p3+g4*p4+g5*p5+g6*p6+g7*p7+g8*p8+g9*p9+g10*p10+g11*p11+g12*p12+g13*p13);;

        err = max( max(r1, r2), max(r3,r4));

        if(err > eps)
        {
            h = h * min( max( facmin, fac*pow(eps/err, 1/9.0)), facmax);

            t = t_prev;
            x = x_prev;
            y = y_prev;
            px = px_prev;
			py = py_prev;
            n--;
        }
        else
        {
            x_prev = x;
            y_prev = y;
            px_prev = px;
            py_prev = py;
            t_prev = t;
            t+=h;
        }

		if(t + h == end)
        {
            cntr=0;
        }
        if((t + h > end)&&(cntr==1))
        {
            h = end - t;
            cntr = 0;
        }
    }
	return;
}

void runge_numbers(double *b)
{
    FILE *number;
    number = fopen("rungenumbers.txt" ,"w");

    for (i = 0; i < 5 ; i++)
    {
        alpha = vec[i];

        for(j=0;j<3;j++)
        {
            xo = 0.0;
            t = 0.0;
            end = M_PI/2;
            ep = 1.e-12;
            eps = (1.e-8)/pow(10,2*j);
            precision = eps;

            h=1;

            if (alpha == 0)
            {
                pxo = 0;
                pyo = 4/ (M_PI - 4);
            }

            yo = -pyo;

            n = 0;
            m = 0;
            q = 0;

            runge(xo, yo, pxo, pyo);

            px_n = px;
            py_n = py;
            x_n = x;
            y_n = y;

           while ((fabs(x-1.0) > precision) || (fabs(py) > precision))
            {

			runge(xo, yo, pxo+ep, pyo);
			a[0][0] = (px - px_n)/ep;
			a[1][0] = (py - py_n)/ep;

			runge(xo, yo, pxo, pyo+ep);
			a[0][1] = (px - px_n)/ep;
			a[1][1] = (py - py_n)/ep;

			a[0][2] = -(x_n-1.0) + a[0][0] * pxo + a[0][1] * pyo;
			a[1][2] = -(py_n) + a[1][0] * pxo + a[1][1] * pyo;

			gauss(a);
			pxo = a[0][2];
            pyo = a[1][2];

			runge(xo, yo, pxo, pyo);

			px_n = px;
			py_n = py;
			x_n = x;
			y_n = y;
            m++;
            }
            end = M_PI/2;
            runge4(xo, yo, pxo, pyo);
            b[j] = x;
            b[j+3] = y;
            b[j+6]= px;
            b[j+9]= py;
        }

        Rx = fabs((b[0] - b[1])/(b[1] - b[2]));
        Ry = fabs((b[3] - b[4])/(b[4] - b[5]));
        Rpx = fabs((b[6] - b[7])/(b[7] - b[8]));
        Rpy = fabs((b[9] - b[10])/(b[10] - b[11]));

        fprintf(number,"\n\nalpha = %5.20lf\n", alpha);
        fprintf(number,"\n & %5.15lf & %5.15lf & %5.15lf & %5.15lf \n", Rx, Ry, Rpx, Rpy);

    }
    printf("\nRunge numbers in 'rungenumbers.txt'\n");
    return;
}
