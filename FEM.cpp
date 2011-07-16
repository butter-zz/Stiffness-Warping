#include <GLUT/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define WIDTH 200
#define HEIGHT 200

typedef struct xyz
{
   long double x;
   long double y;
   long double z;
} point;

typedef struct simentai
{
    point vertex[4];
   long double q[12];             //変位ベクトル行列
   long double invC[4][4];         //形状関数の情報を持った行列
   long double F[12];             //弾性力
   long double V;                 //体積
   long double poison;            //ポアソン比
   long double young;             //ヤング率
   long double lambda;            //ラメの定数(λ)
   long double myu;               //ラメの定数(μ)
   long double m;                 //質点の重さ
    point  velocity[4];   //点ごとの速度
} simentai;

point rotate[2];
simentai object;

void calcInverse(long double A[][4],long double B[][4])
{
    //まずA[i][i]の最大値を探し、最大値がある行とi行を入れ替え1/A[i][i]をして1にし、ほかの行が0になるようにがんばる
    //init inverse matrix
    for(int i = 0;i < 4;i++)
	for(int t = 0;t < 4;t++)
	    B[i][t] = (i == t)?1.0:0.0;

    int check[4] = {0,0,0,0};
    for(int i = 0;i < 4;i++)
	{
	    //search max
	  long double max = 0;
	    for(int t = i;t < 4;t++)
		{
		    if(max < fabs(A[t][i])) 
			{
			    check[i] = t;
			    max = fabs(A[t][i]);
			}
		}
	    if(max == 0) printf("ゼロだよー\n");
	    //change row
	    if(i != check[i]){
		for(int t = 0;t < 4;t++)
		    {
			long double dummy;
			dummy = A[i][t];
			A[i][t] = A[check[i]][t];
			A[check[i]][t] = dummy;
			dummy = B[i][t];
			B[i][t] = B[check[i]][t];
			B[check[i]][t] = dummy;
		    }
	    }
	    //taikaku to 1  
	   long double pivot = 1 / A[i][i];
	    for(int t = 0;t < 4;t++)
		{
		    A[i][t] *= pivot;
		    B[i][t] *= pivot;
		}
	    //other to 0   
	    for(int t = 0;t < 4;t++)
		{
		    if(i != t)
			{
			long double dummy = A[t][i];
			for(int s = 0;s < 4;s++)
			    {
				A[t][s] -= A[i][s] * dummy;
				B[t][s] -= B[i][s] * dummy;
			    }
			}
	    }
	}
}

void calcDisplaceFunc() 
//変位を関数表記した場合の係数を求める
//また、これにより形状関数も求めることが出来る
{
    long double C[4][4] = {
	{1.0,object.vertex[0].x,object.vertex[0].y,object.vertex[0].z},  
	{1.0,object.vertex[1].x,object.vertex[1].y,object.vertex[1].z},  
	{1.0,object.vertex[2].x,object.vertex[2].y,object.vertex[2].z},  
	{1.0,object.vertex[3].x,object.vertex[3].y,object.vertex[3].z}
    };
    calcInverse(C,object.invC);
}

void calcElasticForce()
{
    long double alphaTq = 
	object.invC[1][0] * object.q[0] + object.invC[1][1] * object.q[3] +
	object.invC[1][2] * object.q[6] + object.invC[1][3] * object.q[9];
    long double betaTq = 
	object.invC[2][0] * object.q[1] + object.invC[2][1] * object.q[4] +
	object.invC[2][2] * object.q[7] + object.invC[2][3] * object.q[10];
    long double gammaTq = 
	object.invC[3][0] * object.q[2] + object.invC[3][1] * object.q[5] +
	object.invC[3][2] * object.q[8] + object.invC[3][3] * object.q[11];
    long double xiTq = 
	object.invC[3][0] * object.q[1] + object.invC[2][0] * object.q[2] + 
	object.invC[3][1] * object.q[4] + object.invC[2][1] * object.q[5] + 
	object.invC[3][2] * object.q[7] + object.invC[2][2] * object.q[8] + 
	object.invC[3][3] * object.q[10] + object.invC[2][3] * object.q[11];
    long double roTq = 
	object.invC[3][0] * object.q[0] + object.invC[1][0] * object.q[2] + 
	object.invC[3][1] * object.q[3] + object.invC[1][1] * object.q[5] + 
	object.invC[3][2] * object.q[6] + object.invC[1][2] * object.q[8] + 
	object.invC[3][3] * object.q[9] + object.invC[1][3] * object.q[11];
    long double faiTq = 
	object.invC[2][0] * object.q[0] + object.invC[1][0] * object.q[1] + 
	object.invC[2][1] * object.q[3] + object.invC[1][1] * object.q[4] + 
	object.invC[2][2] * object.q[6] + object.invC[1][2] * object.q[7] + 
	object.invC[2][3] * object.q[9] + object.invC[1][3] * object.q[10];

    long double lVabgTq = object.lambda * object.V * (alphaTq + betaTq + gammaTq);
    long double mV = object.myu * object.V;
    long double alphaTq_2 = 2 * alphaTq;
    long double betaTq_2 = 2 * betaTq;
    long double gammaTq_2 = 2 * gammaTq;

    long double U_lambda[12];
    long double U_myu[12];
    
    for(int i = 0;i < 12;i += 3)
	{
	    U_lambda[i] = lVabgTq * object.invC[1][i/3];
	    
	    U_lambda[i+1] = lVabgTq * object.invC[2][i/3];
	    
	    U_lambda[i+2] = lVabgTq * object.invC[3][i/3];
	 
	    U_myu[i] = mV * (2 * alphaTq * object.invC[1][i/3] + roTq * object.invC[3][i/3] + faiTq * object.invC[2][i/3]);
	    
	    U_myu[i+1] = mV * (2 * betaTq * object.invC[2][i/3] + xiTq * object.invC[3][i/3] + faiTq * object.invC[1][i/3]); 
	    
	    U_myu[i+2] = mV * (2 * gammaTq * object.invC[3][i/3] + xiTq * object.invC[2][i/3] + roTq * object.invC[1][i/3]);
	}

    for(int i = 0;i < 12;i++)
	object.F[i] = -(U_lambda[i] + U_myu[i]);
   
    static int hehe;
    static double counter;
    if(hehe % 100 == 0){
    FILE *fp = fopen("F.txt","a");
    fprintf(fp,"%lf ",counter);
    for(int i = 0;i < 12;i++){
	object.F[i] = -(U_lambda[i] + U_myu[i]);
	fprintf(fp,"%Lf ",object.F[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
    }
    counter += 0.001;
    hehe++;
   }  

void runge()
{
    long double k1[2][4][3];
    long double k2[2][4][3];
    long double k3[2][4][3];
    long double k4[2][4][3];
    long double dt = 0.001;
    
    for(int i = 0;i < 4;i++)
        {
	    k1[0][i][0] = dt * object.velocity[i].x;
	    k1[0][i][1] = dt * object.velocity[i].y;
	    k1[0][i][2] = dt * object.velocity[i].z;
	    
	    k1[1][i][0] = dt * (object.F[i*3]/object.m);
	    k1[1][i][1] = dt * (object.F[i*3+1]/object.m);
	    k1[1][i][2] = dt * (object.F[i*3+2]/object.m);
	    
	    k2[0][i][0] = dt * (object.velocity[i].x + k1[1][i][0] / 2.0);
	    k2[0][i][1] = dt * (object.velocity[i].y + k1[1][i][1] / 2.0);
	    k2[0][i][2] = dt * (object.velocity[i].z + k1[1][i][2] / 2.0);
	    
	    k2[1][i][0] = dt * (object.F[i*3]/object.m + k1[0][i][0] / 2.0);
	    k2[1][i][1] = dt * (object.F[i*3+1]/object.m + k1[0][i][1] / 2.0);
	    k2[1][i][2] = dt * (object.F[i*3+2]/object.m + k1[0][i][2] / 2.0);
	    
	    
	    k3[0][i][0] = dt * (object.velocity[i].x + k2[1][i][0] / 2.0);
	    k3[0][i][1] = dt * (object.velocity[i].y + k2[1][i][1] / 2.0);
	    k3[0][i][2] = dt * (object.velocity[i].z + k2[1][i][2] / 2.0);
	    
	    k3[1][i][0] = dt * (object.F[i*3]/object.m +  k2[0][i][0] / 2.0);
	    k3[1][i][1] = dt * (object.F[i*3+1]/object.m +  k2[0][i][1] / 2.0);
	    k3[1][i][2] = dt * (object.F[i*3+2]/object.m +  k2[0][i][2] / 2.0);
	    
	    k4[0][i][0] = dt * (object.velocity[i].x + k3[1][i][0]);
	    k4[0][i][1] = dt * (object.velocity[i].y + k3[1][i][1]);
	    k4[0][i][2] = dt * (object.velocity[i].z + k3[1][i][2]);
	    
	    k4[1][i][0] = dt * (object.F[i*3]/object.m + k3[0][i][0]);
	    k4[1][i][1] = dt * (object.F[i*3+1]/object.m + k3[0][i][1]);
	    k4[1][i][2] = dt * (object.F[i*3+2]/object.m + k3[0][i][2]);
	    
	    object.q[i*3] += (k1[0][i][0] + 2.0 * k2[0][i][0] + 2.0 * k3[0][i][0] + k4[0][i][0]) / 6.0;
	    object.q[i*3+1] += (k1[0][i][1] + 2.0 * k2[0][i][1] + 2.0 * k3[0][i][1] + k4[0][i][1]) / 6.0;
	    object.q[i*3+2] += (k1[0][i][2] + 2.0 * k2[0][i][2] + 2.0 * k3[0][i][2] + k4[0][i][2]) / 6.0;
	    object.velocity[i].x += (k1[1][i][0] + 2.0 * k2[1][i][0] + 2.0 * k3[1][i][0] + k4[1][i][0]) / 6.0;
	    object.velocity[i].y += (k1[1][i][1] + 2.0 * k2[1][i][1] + 2.0 * k3[1][i][1] + k4[1][i][1]) / 6.0;
	    object.velocity[i].z += (k1[1][i][2] + 2.0 * k2[1][i][2] + 2.0 * k3[1][i][2] + k4[1][i][2]) / 6.0;
	    
        }   
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);

    glLoadIdentity();
    glRotated(rotate[1].y,1,0,0);
    glRotated(rotate[1].x,0,1,0);
    glRotated(0,0,0,1);

    glPointSize(7);
    glBegin(GL_POINTS);
    for(int i = 0;i < 4;i++)
	{
	    if(i == 0) glColor3d(1,0,0);
	    else if(i == 1) glColor3d(0,1,0);
	    else if(i == 2) glColor3d(0,0,1);
	    else if(i == 3) glColor3d(1,0,1);
	    glVertex3d(object.vertex[i].x + object.q[i*3],object.vertex[i].y + object.q[i*3+1],object.vertex[i].z + object.q[i*3+2]);
	}
    glEnd();

    glLineWidth(3);
    glBegin(GL_LINES);
    glColor3d(1,1,0);
    for(int i = 0;i < 4;i++)
	{
	    for(int t = i+1;t < 4;t++)
		{
		    glVertex3d(object.vertex[i].x + object.q[i*3],object.vertex[i].y + object.q[i*3+1],object.vertex[i].z + object.q[i*3+2]);
		    glVertex3d(object.vertex[t].x + object.q[t*3],object.vertex[t].y + object.q[t*3+1],object.vertex[t].z + object.q[t*3+2]);
      		}
	}
    glEnd();
    /*
    static int hehe1;   
    static double counter1;
    if(hehe1 % 100 == 0){
    FILE *fp = fopen("D.txt","a");
    fprintf(fp,"%lf ",counter1);
    for(int i = 0;i < 12;i++){
	fprintf(fp,"%Lf ",object.q[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
    }
    */
    //重心の。。。
/*
    long double g[3] = {0,0,0};
    for(int i = 0;i < 4;i++){
	g[0] += object.vertex[i].x + object.q[i*3];
	g[1] += object.vertex[i].y + object.q[i*3+1];
	g[2] += object.vertex[i].z + object.q[i*3+2];
	//	printf("%Lf %Lf %Lf\n",object.q[i*3],object.q[i*3+1],object.q[i*3+2]);
    }
*/
    //    printf("%Lf %Lf %Lf\n",g[0],g[1],g[2]);
    /*
      if(hehe1 % 100 == 0){
	FILE *fp = fopen("G.txt","a");
	fprintf(fp,"%lf %Lf %Lf %Lf\n",counter1,g[0]/4,g[1]/4,g[2]/4);
	fclose(fp);
	}

    counter1 += 0.001;
    hehe1++;
    */
    calcElasticForce();        
    runge();

    glLineWidth(1);
    glBegin(GL_LINES);
    glColor3d(1,0,0);
    glVertex3d(-100,0,0);
    glVertex3d(100,0,0);
    glColor3d(0,1,0);
    glVertex3d(0,-100,0);
    glVertex3d(0,100,0);
    glColor3d(0,0,1);
    glVertex3d(0,0,-100);
    glVertex3d(0,0,100);
    glEnd();

    glutSwapBuffers();
}

void init()
{
    glClearColor(1,1,1,1);

    object.vertex[0].x = 0;
    object.vertex[0].y = 0; 
    object.vertex[0].z = 0;

    object.vertex[1].x = 1;
    object.vertex[1].y = 0; 
    object.vertex[1].z = 0;

    object.vertex[2].x = 0;
    object.vertex[2].y = 1; 
    object.vertex[2].z = 0;

    object.vertex[3].x = 0;
    object.vertex[3].y = 0; 
    object.vertex[3].z = 1;

    object.q[0] = 0;//0.2;
    object.q[1] = 0;//0.2;
    object.q[2] = 0;//0.2;

    object.q[3] = 0;
    object.q[4] = 0;
    object.q[5] = 0;

    object.q[6] = 0;
    object.q[7] = 0;
    object.q[8] = 0;

    object.q[9] = 0;
    object.q[10] = 0;
    object.q[11] = 0;

    object.poison = 0.4;
    object.young = 3;
    object.lambda = (object.young * object.poison) / ((1 + object.poison) * (1 - 2 * object.poison));
    object.myu = object.young / (2 * (1 + object.poison));
    object.m = 1;

    for(int i = 0;i < 4;i++)
	{
	    object.velocity[i].x = 0;
	    object.velocity[i].y = 0;
	    object.velocity[i].z = 0;
	}

    object.V = fabs(
	(object.vertex[1].x * object.vertex[2].y * object.vertex[3].z + 
	 object.vertex[3].x * object.vertex[1].y * object.vertex[2].z + 
	 object.vertex[2].x * object.vertex[3].y * object.vertex[1].z -
	 object.vertex[3].x * object.vertex[2].y * object.vertex[1].z -
	 object.vertex[2].x * object.vertex[1].y * object.vertex[3].z -
	 object.vertex[1].x * object.vertex[3].y * object.vertex[2].z -
	 object.vertex[0].x * object.vertex[2].y * object.vertex[3].z -
	 object.vertex[3].x * object.vertex[0].y * object.vertex[2].z -
	 object.vertex[2].x * object.vertex[3].y * object.vertex[0].z +
	 object.vertex[3].x * object.vertex[2].y * object.vertex[0].z +
	 object.vertex[2].x * object.vertex[0].y * object.vertex[3].z +
	 object.vertex[0].x * object.vertex[3].y * object.vertex[2].z +
	 object.vertex[0].x * object.vertex[1].y * object.vertex[3].z +
	 object.vertex[1].x * object.vertex[3].y * object.vertex[0].z +
	 object.vertex[3].x * object.vertex[0].y * object.vertex[1].z -
	 object.vertex[3].x * object.vertex[1].y * object.vertex[0].z -
	 object.vertex[1].x * object.vertex[0].y * object.vertex[3].z -
	 object.vertex[0].x * object.vertex[3].y * object.vertex[1].z -
	 object.vertex[0].x * object.vertex[1].y * object.vertex[2].z -
	 object.vertex[2].x * object.vertex[0].y * object.vertex[1].z -
	 object.vertex[1].x * object.vertex[2].y * object.vertex[0].z +
	 object.vertex[2].x * object.vertex[1].y * object.vertex[0].z +
	 object.vertex[1].x * object.vertex[0].y * object.vertex[2].z +
	 object.vertex[0].x * object.vertex[2].y * object.vertex[1].z) / 6.0);

    calcDisplaceFunc();

    for(int i = 0;i < 2;i++)
	{
	    rotate[i].x = 0;
	    rotate[i].y = 0;
	}
}
  
void mouse(int button,int state,int x,int y)
{
    switch(button)
	{
	case GLUT_LEFT_BUTTON:
	    if(state == GLUT_UP);
	    else{
		rotate[0].x = x;
		rotate[0].y = y;
	    }                                          
	    break;
	case GLUT_RIGHT_BUTTON:
	    break;
	case GLUT_MIDDLE_BUTTON:
	    break;
	}
}

void idle(void)
{
    glutPostRedisplay();
}

void motion(int x,int y)
{
    rotate[1].x = x - rotate[0].x;
    rotate[1].y = y - rotate[0].y;
    glutPostRedisplay();
}
void keyboard(unsigned char key,int x,int y)
{
    switch(key)
	{
	case 'p': 
	    object.q[0] = 1;
	    object.q[1] = 1;
	    object.q[2] = 1;
	    break;
	case 'o':
	    object.q[3] = 1;
	    object.q[4] = 1;
	    object.q[5] = -1;
	    break;
	case 'r': 
	    object.vertex[0].x = 0;
	    object.vertex[0].y = 0; 
	    object.vertex[0].z = 0;

	    object.vertex[1].x = 3;
	    object.vertex[1].y = 0; 
	    object.vertex[1].z = 0;

	    object.vertex[2].x = 0;
	    object.vertex[2].y = 3; 
	    object.vertex[2].z = 0;

	    object.vertex[3].x = 0;
	    object.vertex[3].y = 0; 
	    object.vertex[3].z = 3;

	    for(int i = 0;i < 4;i++){
		object.q[i*3  ] = 0;
		object.q[i*3+1] = 0;
		object.q[i*3+2] = 0;
		object.F[i*3  ] = 0;
		object.F[i*3+1] = 0;
		object.F[i*3+2] = 0;
	    }
	    break;
	case 'c':
	    object.vertex[0].x = 0;
	    object.vertex[0].y = 0; 
	    object.vertex[0].z = 0;

	    object.vertex[1].x = 0;
	    object.vertex[1].y = 0; 
	    object.vertex[1].z = 3;

	    object.vertex[2].x = 0;
	    object.vertex[2].y = 3; 
	    object.vertex[2].z = 3;

	    object.vertex[3].x = 3;
	    object.vertex[3].y = 0; 
	    object.vertex[3].z = 3;

	    break;
	case 'q':
	case 'Q':
	case '\033': exit(0); break;
	default: break;
	}
}

void resize(int w,int h)                                              
{
    glViewport(0,0,w,h); //ウィンドウ全体をビューポートに        
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    //スクリーン上の表示領域をビューポートに比例させる                    
    glOrtho(-5, 5, -5,5, -5,5);
    //    glOrtho( - 0.5, (GLdouble)w -0.5, (GLdouble)h -0.5, -0.5,-1,1); 
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char *argv[])
{
    glutInitWindowPosition(700,500);
    glutInitWindowSize(400,320);
    glutInit(&argc,argv);
    glutCreateWindow(argv[0]);
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE);
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutReshapeFunc(resize);
    init();
    glutMainLoop();
    return 0;
}


