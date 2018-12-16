/**********************************************************************
 * Try the following command to compile and run this program:
 * % ./execute.sh
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <GL/glut.h>

static double interval = 0.001;
static double cur_t = 0;
static double tot_t = 10;

static double REF_X = 0.5;
static double REF_Y = 0.3;

/* stiffness coefficient */
static double K_X = 10;
static double K_Y = 3;
/* joint damping coefficient */
static double D = 0.1;
/* damping coefficient */
static double C_X = 2;
static double C_Y = 2.2;

static double ERROR_P = 1e-3;
static double ERROR_V = 1e-2;

FILE *endpoint_fp;
FILE *joint_fp;

/**********************************************************************
 arm section
 **********************************************************************/

double linklength[] = {
  0.05,
  0.10,
  0.25,
  0.25,
  0.15,
};

double linkmass[] = {
  0.0,
  0.5,
  1.0,
  1.0,
  0.5,
};

double linkcomx[] = {
  0.0,
  0.05,
  0.10,
  0.10,
  0.05,
};

double linkcomy[] = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
};

double linkinertia[] = {
  0.0,
  0.005,
  0.01,
  0.007,
  0.003,
};

/* given values for equation of motion */
double tau[] = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
};

double fx = 0.0;

double fy = 0.0;

double q[] = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
};

double dq[] = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
};

/**********************************************************************
 misc
 **********************************************************************/

/*! \brief generate a number between min and max at random */
double frand(double min, double max)
{
  return min + ( max - min ) * ( (double)rand() / RAND_MAX );
}

/*! \brief calculate max within two values  */
int max_val(int a, int b)
{
  if ( a > b ){
    return a;
  }
  else {
    return b;
  }
}

/**********************************************************************
 vector-matrix algebra section
 **********************************************************************/

/*! \brief copy a vector */
void vec_copy(double src[], double dest[], int n)
{
  memcpy( dest, src, sizeof(double)*n );
}

/*! \brief add two vectors */
void vec_add(double v1[], double v2[], double v[], int n)
{
  for( ; n>0; n-- )
    *v++ = *v1++ + *v2++;
}

/*! \brief subtract a vector from another */
void vec_sub(double v1[], double v2[], double v[], int n)
{
  for( ; n>0; n-- )
    *v++ = *v1++ - *v2++;
}

/*! \brief inner product of two vectors */
double vec_innerprod(double v1[], double v2[], int n)
{
  double val = 0;

  for( ; n>0; n-- )
    val += *v1++ * *v2++;
  return val;
}

/*! \brief generate a vector at random */
void vec_rand(double v[], int n, double min, double max)
{
  for( ; n>0; n-- )
    *v++ = frand( min, max );
}

/*! \brief write a vector to the standard output */
void vec_write(double v[], int n)
{
  int i;

  for( i=0; i<n; i++ )
    printf( " %lf", v[i] );
  printf( "\n" );
}

/*! \brief copy a matrix */
void mat_copy(double src[], double dest[], int nr, int nc)
{
  memcpy( dest, src, sizeof(double)*nr*nc );
}

/*! \brief generate a matrix at random */
void mat_rand(double m[], int nr, int nc, double min, double max)
{
  int n;

  n = nr * nc;
  for( ; n>0; n-- )
    *m++ = frand( min, max );
}

/*! \brief write a matrix to the standard output */
void mat_write(double m[], int nr, int nc)
{
  int i, j;

  for( i=0; i<nr; i++, m+=nc ){
    for( j=0; j<nc; j++ )
      printf( " %lf", *(m+j) );
    printf( "\n" );
  }
}

/*! \brief multiply a matrix to a vector from the left side */
void mat_vec_mul(double m[], double v1[], double v[], int n)
{
  int i;

  for( i=0; i<n; i++, m+=n )
    *v++ = vec_innerprod( m, v1, n );
}

/*! \brief solve a linear equation */
void le_solve(double a[], double b[], int n)
{
  int i, j, k;

  for ( k=1; k<=n; k++ ){
    for ( i=k+1; i<=n; i++ ){
      b[i-1] -= a[(i-1)*n+k-1] / a[(k-1)*n+k-1] * b[k-1];
      for ( j=n; j>=k; j-- ){
        a[(i-1)*n+j-1] -= a[(i-1)*n+k-1] / a[(k-1)*n+k-1] * a[(k-1)*n+j-1];
      }
    }
  }

  for ( k=n; k>=1; k-- ){
    for ( i=k-1; i>=1; i-- ){
      b[i-1] -= a[(i-1)*n+k-1] / a[(k-1)*n+k-1] * b[k-1];
      a[(i-1)*n+k-1] = 0;
    }
    b[k-1] /= a[(k-1)*n+k-1];
    a[(k-1)*n+k-1] = 1;
  }
}

/**********************************************************************
 arm model section
 **********************************************************************/

/*! \brief a link with a revolute joint */
typedef struct{
  double l;     /*!< link length */
  double m;     /*!< mass */
  double rgx;   /*!< longitudinal distance of the center of mass */
  double rgy;   /*!< leftward offset of the center of mass */
  double i;     /*!< moment of inertia about the center of mass */

  double q;     /*!< link joint angle */
  double dq;    /*!< link joint velocity */
  double ddq;   /*!< link joint acceleration */

  double _theta;   /* accumulated joint angle */
  double _dtheta;  /* accumulated joint velocity */
  double _ddtheta; /* accumulated joint acceleration */
  double _x;       /* x-position of joint */
  double _y;       /* y-position of joint */
  double _vx;      /* x-velocity of joint */
  double _vy;      /* y-velocity of joint */
  double _ax;      /* x-acceleration of joint */
  double _ay;      /* y-acceleration of joint */
  double _xg;      /* x-position of the center of mass */
  double _yg;      /* y-position of the center of mass */
  double _agx;     /* x-acceleration of the center of mass */
  double _agy;     /* y-acceleration of the center of mass */
  double _jx;      /* basis of Jacobian in x */
  double _jy;      /* basis of Jacobian in y */
  double _hx;      /* basis of centripetal / Coriolis acc. in x */
  double _hy;      /* basis of centripetal / Coriolis acc. in y */

  double _fx;      /* x-constraint force at joint */
  double _fy;      /* y-constraint force at joint */
  double tau;      /* actuation torque at joint */
  double fex;      /* x-external force */
  double fey;      /* y-external force */
  double xe;       /* x-position of point of action of external force */
  double ye;       /* y-position of point of action of external force */

  double _lc;      /* l cos(theta) */
  double _ls;      /* l sin(theta) */
} link_t;

/*! \brief create a link. */
void link_create(link_t *link, double l, double m, double rgx, double rgy, double i)
{
  link->l = l;
  link->m = m;
  link->rgx = rgx;
  link->rgy = rgy;
  link->i = i;
  link->q = link->dq = link->ddq = 0;

  link->_theta = link->_dtheta = link->_ddtheta = 0;
  link->_x = l;
  link->_y = 0;
  link->_vx = link->_vy = 0;
  link->_ax = link->_ay = 0;
  link->_xg = link->_yg = 0;
  link->_agx = link->_agy = 0;
  link->_jx = link->_jy = 0;
  link->_hx = link->_hy = 0;

  link->_fx = link->_fy = 0;
  link->tau = 0;
  link->fex = link->fey = 0;
  link->xe = link->ye = 0;
}

/*! \brief skeleton */
typedef struct skeleton_t {
  int num;      /*!< number of links */
  link_t *link; /*!< root link */

  double *H;    /* inertia matrix */
  double *b;    /* bias force vector */
  double *J;    /* transpose of Jacobian matrix */
  double *f;    /* external force */
} skeleton_t;

/* brief calculate h_xk */
double h_xk(skeleton_t *skeleton, int k, int j)
{
  skeleton_t *s = skeleton;
  double h_xk;

  if( k < j ){
    h_xk = 0;
  }
  else if( k == j ){
    h_xk = s->link[j].rgx * cos( ( &s->link[j] )->_theta )
         - s->link[j].rgy * sin( ( &s->link[j] )->_theta );
  }
  else {
    h_xk = s->link[j].l * cos( ( &s->link[j] )->_theta );
  }

  return h_xk;
}

/* brief calculate h_yk */
double h_yk(skeleton_t *skeleton, int k, int j)
{
  skeleton_t *s = skeleton;
  double h_yk;

  if( k < j ){
    h_yk = 0;
  }
  else if( k == j ){
    h_yk = s->link[j].rgx * sin( ( &s->link[j] )->_theta )
         - s->link[j].rgy * cos( ( &s->link[j] )->_theta );
  }
  else {
    h_yk = s->link[j].l * sin( ( &s->link[j] )->_theta );
  }

  return h_yk;
}

/*! \brief update position, velocity and acceleration of a skeleton */
void skeleton_update_state(skeleton_t *skeleton)
{
  int i;
  link_t *l, *lp;

  lp = &skeleton->link[0];
  for( i=1; i<skeleton->num; i++, lp=l ){
    l = &skeleton->link[i];
    /* position */
    l->_theta = lp->_theta + l->q;
    l->_x = lp->_x + ( l->_lc = l->l * cos( l->_theta ) );
    l->_y = lp->_y + ( l->_ls = l->l * sin( l->_theta ) );
    /* position of center of mass */
    l->_xg = lp->_x + l->rgx * cos(l->_theta) - l->rgy * sin(l->_theta);
    l->_yg = lp->_y + l->rgx * sin(l->_theta) + l->rgy * cos(l->_theta);
    /* velocity */
    l->_dtheta = lp->_dtheta + l->dq;
    l->_vx = lp->_vx - l->_dtheta * l->_ls;
    l->_vy = lp->_vy + l->_dtheta * l->_lc;
    /* acceleration */
    l->_ddtheta = lp->_ddtheta + l->ddq;
    l->_ax = lp->_ax
      - l->_dtheta * l->_dtheta * l->_lc - l->_ddtheta * l->_ls;
    l->_ay = lp->_ay
      - l->_dtheta * l->_dtheta * l->_ls + l->_ddtheta * l->_lc;
    /* acceleration of center of mass */
    l->_agx = lp->_ax - l->_ddtheta * (l->rgx * sin(l->_theta) + l->rgy * cos(l->_theta)) - pow(l->_dtheta, 2) * (l->rgx * cos(l->_theta) - l->rgy * sin(l->_theta));
    l->_agy = lp->_ay + l->_ddtheta * (l->rgx * cos(l->_theta) - l->rgy * sin(l->_theta)) - pow(l->_dtheta, 2) * (l->rgx * sin(l->_theta) - l->rgy * cos(l->_theta));
  }
}

/*! \brief update inertia matrix, system bias force and Jacobian matrix of a skeleton */
void skeleton_update_fdstate(skeleton_t *skeleton)
{
  int i, j, k;
  int n = skeleton->num;
  skeleton_t *s = skeleton;

  /* compute H */
  for( i=1; i<n; i++ ){
    for( j=1; j<n; j++ ){
      s->H[(i-1)*(n-1)+j-1] = 0;
      for( k=max_val( i, j ); k<n; k++ ){
        s->H[(i-1)*(n-1)+j-1] += s->link[k].m * ( ( s->link[k]._xg - s->link[i-1]._x ) * ( s->link[k]._xg - s->link[j-1]._x )
                                                + ( s->link[k]._yg - s->link[i-1]._y ) * ( s->link[k]._yg - s->link[j-1]._y ) ) + s->link[k].i;
      }
    }
  }

  /* compute b */
  double bb;
  for( i=1; i<n; i++ ){
    s->b[i-1] = 0;
    for( k=i; k<n; k++ ){
      bb = 0;
      for( j=1; j<n; j++ ){
        bb += ( ( s->link[k]._yg - s->link[i-1]._y ) * h_xk( s, k, j )
              - ( s->link[k]._xg - s->link[i-1]._x ) * h_yk( s, k, j ) )
              * ( &s->link[j] )->_dtheta * ( &s->link[j] )->_dtheta;
      }
      s->b[i-1] +=  s->link[k].m * bb;
    }
  }

  /* compute J_E^T */
  for( i=1; i<n; i++ ){
    s->J[(i-1)*2]   = - ( s->link[n-1].ye - s->link[i-1]._y );
    s->J[(i-1)*2+1] =     s->link[n-1].xe - s->link[i-1]._x  ;
  }
}

/**********************************************************************
 inverse dynamics section
 **********************************************************************/

/*! \brief update joint forces of a skeleton */
void skeleton_update_force(skeleton_t *skeleton)
{
  int i;
  double fgx, fgy; /* net force on a link */
  link_t *l, *lc, *lp; /* link, child link, parent link */

  lc = NULL;
  for( i=skeleton->num-1; i>0; i--, lc=l ){
    l = &skeleton->link[i];
    lp = &skeleton->link[i-1];
    fgx = l->m * l->_agx;
    fgy = l->m * l->_agy;
    l->_fx = fgx - l->fex;
    l->_fy = fgy - l->fey;
    l->tau = l->i * l->_ddtheta - (l->xe - lp->_x) * l->fey + (l->ye - lp->_y) * l->fex + (l->_xg - lp->_x) * fgy - (l->_yg - lp->_y) * fgx;
    if( lc != NULL ){
      l->_fx += lc->_fx;
      l->_fy += lc->_fy;
      l->tau += lc->tau + (l->_x - lp->_x) * lc->_fy - (l->_y - lp->_y) * lc->_fx;
    }
  }

  /* compute f_E */
  skeleton->f[0] = skeleton->link[skeleton->num-1].fex;
  skeleton->f[1] = skeleton->link[skeleton->num-1].fey;
}

/**********************************************************************
 forward dynamics section
 **********************************************************************/

/* \brief calculate Jaconian matrix multiply external force */
void compute_Jf(skeleton_t *skeleton, double Jf[]){
  int i;
  for( i=1; i<skeleton->num; i++ ){
    Jf[i-1] = skeleton->J[(i-1)*2] * skeleton->f[0] + skeleton->J[(i-1)*2+1] * skeleton->f[1];
  }
}

/*! \brief calculate Jf and solve forward dynamics of a skeleton */
void skeleton_fd(skeleton_t *skeleton, double tau[])
{
  int n = skeleton->num;
  double *Jf;
  double *v1;
  double *v2;

  if( ( Jf = malloc( sizeof(double)*(n-1) ) ) == NULL ){
    fprintf( stderr, "cannot allocate memory for links.\n" );
    return;
  }
  if( ( v1 = malloc( sizeof(double)*(n-1) ) ) == NULL ){
    fprintf( stderr, "cannot allocate memory for links.\n" );
    return;
  }
  if( ( v2 = malloc( sizeof(double)*(n-1) ) ) == NULL ){
    fprintf( stderr, "cannot allocate memory for links.\n" );
    return;
  }

  compute_Jf( skeleton, Jf );

  /* compute v1 = b - J_E^T * f_E */
  vec_sub( skeleton->b, Jf, v1, n-1 );

  /* compute v2 = tau - b + J_E^T * f_E */
  vec_sub( &tau[1], v1, v2, skeleton->num-1 );

  /* compute v2 = ddq */
  le_solve( skeleton->H, v2, skeleton->num-1 );

  printf( "H:\n" );
  mat_write( skeleton->H, skeleton->num-1, skeleton->num-1 );
  printf( "f:\n" );
  vec_write( skeleton->f, 2 );
  printf( "b:\n" );
  vec_write( skeleton->b, skeleton->num-1 );
  printf( "v1:\n" );
  vec_write( v1, skeleton->num-1 );
  printf( "tau:\n" );
  vec_write( &tau[1], skeleton->num-1 );
  printf( "v2:\n" );
  vec_write( v2, skeleton->num-1 );

  int i;
  for ( i=1; i<skeleton->num; i++ ){
    skeleton->link[i].ddq = v2[i-1];
  }

  free( Jf );
  free( v1 );
  free( v2 );
}

/*! \brief create a skeleton model from an array of link lengths */
skeleton_t *skeleton_create(skeleton_t *skeleton, double l[], double m[], double rgx[], double rgy[], double in[], int n)
{
  if( ( skeleton->link = malloc( sizeof(link_t)*n ) ) == NULL ){
    fprintf( stderr, "cannot allocate memory for links.\n" );
    return NULL;
  }
  if( ( skeleton->H = malloc( sizeof(double)*(n-1)*(n-1) ) ) == NULL ){
    fprintf( stderr, "cannot allocate memory for links.\n" );
    return NULL;
  }
  if( ( skeleton->b = malloc( sizeof(double)*(n-1) ) ) == NULL ){
    fprintf( stderr, "cannot allocate memory for links.\n" );
    return NULL;
  }
  if( ( skeleton->J = malloc( sizeof(double)*(n-1)*2 ) ) == NULL ){
    fprintf( stderr, "cannot allocate memory for links.\n" );
    return NULL;
  }
  if( ( skeleton->f = malloc( sizeof(double)*2 ) ) == NULL ){
    fprintf( stderr, "cannot allocate memory for links.\n" );
    return NULL;
  }

  int i;
  skeleton->num = n;
  for( i=0; i<n; i++ ){
    link_create( &skeleton->link[i], l[i], m[i], rgx[i], rgy[i], in[i] );
  }

  skeleton_update_state( skeleton );
  skeleton_update_fdstate( skeleton );
  skeleton_update_force( skeleton );

  return skeleton;
}

/**********************************************************************
 Controller section
 **********************************************************************/

void control_simple( skeleton_t *skeleton, double ref_x, double ref_y, double q[], double dq[] )
{
  double Fx, Fy;

  Fx = K_X * ( ref_x - skeleton->link[skeleton->num-1]._x ) - C_X * skeleton->link[skeleton->num-1]._vx;
  Fy = K_Y * ( ref_y - skeleton->link[skeleton->num-1]._y ) - C_Y * skeleton->link[skeleton->num-1]._vy;

  /* printf( "Fx:%f, Fy:%f\n", Fx, Fy ); */
  /* getchar(); */

  int i;

  /* input q, dq, fx and fy */
  for( i=1; i<skeleton->num; i++ ){
    skeleton->link[i].q   = q[i];
    skeleton->link[i].dq  = dq[i];
    skeleton->link[i].ddq = 0;
  }

  skeleton->link[skeleton->num-1].xe = skeleton->link[skeleton->num-1]._x;
  skeleton->link[skeleton->num-1].ye = skeleton->link[skeleton->num-1]._y;
  skeleton->link[skeleton->num-1].fex = -Fx;
  skeleton->link[skeleton->num-1].fey = -Fy;

  skeleton_update_state( skeleton );
  skeleton_update_force( skeleton );

  /* output q, dq, tau */
  printf( "Fx: %f Fy: %f\n", Fx, Fy );
  printf( "control_simple_tau:" );
  for( i=1; i<skeleton->num; i++ ){
    tau[i] = skeleton->link[i].tau - D * skeleton->link[i].dq;
    printf( "%f ", tau[i] );
  }
  printf( "\n" );
}

/**********************************************************************
 Approximation section
 **********************************************************************/

void euler( skeleton_t *skeleton, double tau[], double q[], double dq[], double ddq[] )
{
  int i;

  for( i=0; i<skeleton->num; i++ ){
     q[i] += interval *  dq[i];
    dq[i] += interval * ddq[i];
  }
}

void r_slope( skeleton_t *skeleton, double tau[], double q[], double dq[], double ddq[] )
{
  int i;
  int n = skeleton->num;
  double k1[n], k2[n], k3[n], k4[n];
  double k5[n], k6[n], k7[n], k8[n];
  double tmpdq[n], tmpq[n];

  /* compute k1, k5 */
  vec_copy( ddq, k1, n );
  vec_copy(  dq, k5, n );

  /* compute k2, k6 */
  vec_copy( dq, tmpdq, n );
  vec_copy(  q,  tmpq, n );

  for( i=0; i<n; i++ ){
    tmpdq[i] = dq[i] + interval / 2 * k1[i];
     tmpq[i] =  q[i] + interval / 2 * k5[i];
  }

  for( i=0; i<n; i++ ){
    skeleton->link[i].q  =  tmpq[i];
    skeleton->link[i].dq = tmpdq[i];
  }

  skeleton_update_state( skeleton );

  skeleton->link[skeleton->num-1].xe = skeleton->link[skeleton->num-1]._x;
  skeleton->link[skeleton->num-1].ye = skeleton->link[skeleton->num-1]._y;

  skeleton_update_fdstate( skeleton );
  skeleton_update_force( skeleton );

  skeleton_fd( skeleton, tau );

  for( i=0; i<n; i++ ){
    k2[i] = skeleton->link[i].ddq;
    k6[i] = skeleton->link[i].dq;
  }

  /* compute k3, k7 */
  for( i=0; i<n; i++ ){
    tmpdq[i] = dq[i] + interval / 2 * k2[i];
     tmpq[i] =  q[i] + interval / 2 * k6[i];
  }

  for( i=0; i<n; i++ ){
    skeleton->link[i].q  =  tmpq[i];
    skeleton->link[i].dq = tmpdq[i];
  }

  skeleton_update_state( skeleton );

  skeleton->link[skeleton->num-1].xe = skeleton->link[skeleton->num-1]._x;
  skeleton->link[skeleton->num-1].ye = skeleton->link[skeleton->num-1]._y;

  skeleton_update_fdstate( skeleton );
  skeleton_update_force( skeleton );

  skeleton_fd( skeleton, tau );

  for( i=0; i<n; i++ ){
    k3[i] = skeleton->link[i].ddq;
    k7[i] = skeleton->link[i].dq;
  }

  /* compute k4, k8 */
  for( i=0; i<n; i++ ){
    tmpdq[i] = dq[i] + interval * k3[i];
     tmpq[i] =  q[i] + interval * k7[i];
  }

  for( i=0; i<n; i++ ){
    skeleton->link[i].q  =  tmpq[i];
    skeleton->link[i].dq = tmpdq[i];
  }

  skeleton_update_state( skeleton );

  skeleton->link[skeleton->num-1].xe = skeleton->link[skeleton->num-1]._x;
  skeleton->link[skeleton->num-1].ye = skeleton->link[skeleton->num-1]._y;

  skeleton_update_fdstate( skeleton );
  skeleton_update_force( skeleton );

  skeleton_fd( skeleton, tau );

  for( i=0; i<n; i++ ){
    k4[i] = skeleton->link[i].ddq;
    k8[i] = skeleton->link[i].dq;
  }

  /* compute ddq */
  for( i=0; i<n; i++ ){
    ddq[i] = ( k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i] ) / 6;
     dq[i] = ( k5[i] + 2 * k6[i] + 2 * k7[i] + k8[i] ) / 6;
  }
}

void rk4( skeleton_t *skeleton, double tau[], double q[], double dq[], double ddq[] )
{
  int i;
  double tmpdq[skeleton->num];

  vec_copy( dq, tmpdq, skeleton->num );

  r_slope( skeleton, tau, q, tmpdq, ddq );
  for( i=0; i<skeleton->num; i++ ){
     q[i] += interval * tmpdq[i];
    dq[i] += interval * ddq[i];
  }
}

void skeleton_evolve( skeleton_t *skeleton, double tau[], double fx, double fy, double q[], double dq[] )
{
  int i;
  double ddq[skeleton->num];

  skeleton->link[skeleton->num-1].xe = skeleton->link[skeleton->num-1]._x;
  skeleton->link[skeleton->num-1].ye = skeleton->link[skeleton->num-1]._y;
  skeleton->link[skeleton->num-1].fex = fx;
  skeleton->link[skeleton->num-1].fey = fy;

  skeleton_update_state( skeleton );
  skeleton_update_force( skeleton );

  /* calculate forward dynamics to get ddq */
  skeleton_update_fdstate( skeleton );
  skeleton_fd( skeleton, tau );
  printf( "skeleton_evolve_tau:" );
  vec_write( tau, skeleton->num );

  for( i=0; i<skeleton->num; i++ ){
    ddq[i] = skeleton->link[i].ddq;
  }

  skeleton_update_state( skeleton );

  printf( "before integral ");
  printf( "q   " ); vec_write( q, skeleton->num );
  printf( "dq  " ); vec_write( dq, skeleton->num );
  printf( "ddq " ); vec_write( ddq, skeleton->num );
#ifdef EULER
  euler( skeleton, tau, q, dq, ddq );
#else
  rk4( skeleton, tau, q, dq, ddq );
#endif
  printf( "after integral ");
  printf( "q   " ); vec_write( q, skeleton->num );
  printf( "dq  " ); vec_write( dq, skeleton->num );
  printf( "ddq " ); vec_write( ddq, skeleton->num );

  for( i=0; i<skeleton->num; i++ ){
    skeleton->link[i].dq = dq[i];
    skeleton->link[i]. q =  q[i];
    ddq[i] = skeleton->link[i].ddq;
  }

  skeleton_update_state( skeleton );
  skeleton_update_force( skeleton );
}

/* arm */
skeleton_t arm;

/**********************************************************************
 OpenGL section
 **********************************************************************/

/* window name */
char window_name[BUFSIZ];

/* window size */
#define WINDOW_WIDTH  640
#define WINDOW_HEIGHT 480

/* viewport */
#define WINDOW_XMIN (-0.6)
#define WINDOW_XMAX ( 1.0)
#define WINDOW_YMIN (-0.3)
#define WINDOW_YMAX ( 0.9)

/*! \brief draw a skeleton */
void draw_skeleton(skeleton_t *skeleton)
{
  int i;

  glColor3f( 0.0, 0.0, 1.0 );
  glLineWidth( 3.0 );
  glBegin( GL_LINE_STRIP );
  glVertex2d( 0, 0 );
  for( i=0; i<skeleton->num; i++ ){
    glVertex2d( skeleton->link[i]._x, skeleton->link[i]._y );
  }
  glEnd();
  glColor3f( 0.0, 1.0, 0.0 );
  glPointSize( 5.0 );
  glBegin( GL_POINTS );
  for( i=0; i<skeleton->num; i++ ){
    glVertex2d( skeleton->link[i]._x, skeleton->link[i]._y );
  }
  glEnd();
}

/*! \brief capture snapshots of animation by utilizing xwd and ImageMagick. */
void capture(void)
{
  static int count = 0;
  static char xwdcommand[BUFSIZ];
  static char xwdfile[BUFSIZ], imgfile[BUFSIZ];

  if( count % 20 == 0 ){
    sprintf( xwdfile, "cap%05d.xwd", count );
    sprintf( imgfile, "cap%05d.png", count );
    sprintf( xwdcommand, "xwd -silent -name %s -nobdrs -out %s", window_name, xwdfile );
    system( xwdcommand );
    sprintf( xwdcommand, "convert %s %s", xwdfile, imgfile );
    system( xwdcommand );
  }
  count++;
}

/*! \brief draw event callback function. */
void draw_main(void)
{
  glClear( GL_COLOR_BUFFER_BIT );

  /* drawing the target point */
  glColor3f( 1.0, 0.0, 0.0 );
  glPointSize( 5.0 );
  glBegin( GL_POINTS );
    glVertex2d( REF_X, REF_Y );
  glEnd();

  if( cur_t < tot_t &&
      ( fabs(arm.link[arm.num-1]._x - REF_X) > ERROR_P ||
        fabs(arm.link[arm.num-1]._y - REF_Y) > ERROR_P ||
        fabs(arm.link[arm.num-1]._vx)        > ERROR_V ||
        fabs(arm.link[arm.num-1]._vy)        > ERROR_V ) ){
    control_simple( &arm, REF_X, REF_Y, q, dq );

    printf( "************************\n" );
    printf( "fx:%f fy:%f\n", fx, fy );
    printf( "************************\n" );
    skeleton_evolve( &arm, tau, fx, fy, q, dq );

    cur_t += interval;
    printf( "cur_t = %f\n", cur_t );

    printf( "fabs(REF_X - x) = %f, fabs(REF_Y - y) = %f\n", fabs(REF_X - arm.link[arm.num-1]._x), fabs(REF_Y - arm.link[arm.num-1]._y) );
    printf( "vx = %f, vy = %f\n", arm.link[arm.num-1]._vx, arm.link[arm.num-1]._vy );

    fprintf( endpoint_fp, "%f %f %f %f %f %f %f\n", cur_t, arm.link[arm.num-1]._x, arm.link[arm.num-1]._y, arm.link[arm.num-1]._vx, arm.link[arm.num-1]._vy, arm.link[arm.num-1]._ax, arm.link[arm.num-1]._ay );

    int i;

    fprintf( joint_fp, "%f ", cur_t );
    for( i=0; i<arm.num; i++ ){
      fprintf( joint_fp, "%f ", arm.link[i].q );
    }
    for( i=0; i<arm.num; i++ ){
      fprintf( joint_fp, "%f ", arm.link[i].dq );
    }
    for( i=0; i<arm.num; i++ ){
      fprintf( joint_fp, "%f ", arm.link[i].ddq );
    }
    for( i=0; i<arm.num; i++ ){
      fprintf( joint_fp, "%f ", arm.link[i].tau );
    }
    fprintf( joint_fp, "\n" );
  }
  else{
    return;
  }

  draw_skeleton( &arm );
  glutSwapBuffers();
  /* if you want to capture snapshots, uncomment the following line. */
  capture();
}

/*! \brief initialize glut objects. */
void draw_init(int *argc, char *argv[])
{
  glutInit( argc, argv );
  /* 2D OpenGL drawing doesn't use the depth buffer */
  glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
  /* window geometry */
#ifdef EULER
  glutInitWindowPosition( 0, 0 );
#else
  glutInitWindowPosition( 0, WINDOW_HEIGHT );
#endif
  glutInitWindowSize( WINDOW_WIDTH, WINDOW_HEIGHT );
  strcpy( window_name, argv[0] );
  glutCreateWindow( argv[0] );
  glClearColor( 0.9, 0.9, 0.9, 1.0 );
  /* invariant projection matrix */
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluOrtho2D( WINDOW_XMIN, WINDOW_XMAX, WINDOW_YMIN, WINDOW_YMAX );
  glMatrixMode( GL_MODELVIEW );
}

/*! \brief key event callback function. */
void draw_key(unsigned char key, int x, int y)
{ /* exit program by typing ESC or 'q' */
  if( key == '\033' || key == 'q' ) exit( 0 );
}

/*! \brief idle event callback function. */
void draw_idle(void)
{
  glutPostRedisplay();
}

/**********************************************************************
 main section
 **********************************************************************/

int main(int argc, char *argv[])
{
  srand( (int)time( NULL ) );
  if( skeleton_create( &arm, linklength, linkmass, linkcomx, linkcomy, linkinertia, sizeof(linklength)/sizeof(double) ) == NULL )
    exit( EXIT_FAILURE );

  tau[0] = q[0] = 0.0;

  endpoint_fp = fopen( "endpoint.dat", "w" );
  joint_fp = fopen( "joint.dat", "w" );

  draw_init( &argc, argv );
  glutDisplayFunc( draw_main );
  glutKeyboardFunc( draw_key );
  glutIdleFunc( draw_idle );
  glutMainLoop();

  fclose( endpoint_fp );
  fclose( joint_fp );

  return 0;
}
