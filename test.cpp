
#include <GLUT/glut.h>
#include <iostream>
#define _USE_MATH_DEFINES //M_PIを使用するための定義
#include <math.h>
#include "calculation.h"

using namespace std;

#define WINDOW_WIDTH 300   //ウィンドウの幅
#define WINDOW_HEIGHT 300  //ウィンドウの高さ
#define VIEW_WIDTH 200     //描画範囲の幅
#define VIEW_HEIGHT 200    //描画範囲の高さ

//四面体の点
double vertex[4][3] = {
    {0.5,0,0.5},
    {0.5,1,0},
    {1,0,0},
    {0,0,0},
   };

//四面体のエッジ
int edge[][2] = {
    {0,1},
    {0,2},
    {0,3},
    {1,2},
    {1,3},
    {2,3},
   };

//四面体の点ごとの重さ
double weight[] = {
    1,
    1,
    1,
    1};

//回転後の座標
GLdouble rotateVertex[4][3] = {
    {0,0,0},
    {0,0,0}, 
    {0,0,0},
    {0,0,0}
};

double g[4][3]; //マッチング結果の座標
double r_x = 0,r_y = 0,r_z = 0; //角度(degree)    
double move = 0;

void Inverse(double mat[][3],double inv[][3]);

void Shape_Matching()
{
    ////x軸回転                                                         
    double Xrotate[][4] ={
	{1,0,0,0},
	{0,cos(r_x * (M_PI / 180)),-sin(r_x * (M_PI / 180)),0},                       
	{0,sin(r_x * (M_PI / 180)),cos(r_x * (M_PI / 180)),0},
	{0,0,0,1}
    };

    ////y軸回転                                                                     
    double Yrotate[][4] ={
	{cos(r_y * (M_PI / 180)),0,sin(r_y * (M_PI / 180)),0},
	{0,1,0,0},
	{-sin(r_y * (M_PI / 180)),0,cos(r_y * (M_PI / 180)),0},
	{0,0,0,1}
    };
    
    //z軸回転                                                                       
    double Zrotate[][4] ={
	{cos(r_z * (M_PI / 180)),-sin(r_z * (M_PI / 180)),0,0},
	{sin(r_z * (M_PI / 180)),cos(r_z * (M_PI / 180)),0,0},
	{0,0,1,0},
	{0,0,0,1}
    };

    double subver[4][3]; //回転計算に使用
    
    for(int i = 0;i < 4;i++)
        for(int t = 0;t < 3;t++)
            subver[i][t] = vertex[i][t] * 2;

    //回転後の座標を計算
        for(int i = 0;i < 4;i++)                                          
    	for(int t = 0;t < 3;t++)                                        
	    rotateVertex[i][t] = Yrotate[t][0] * subver[i][0] + Yrotate[t][1] * subver[i][1] + Yrotate[t][2] * subver[i][2];
    
    for(int i = 0;i < 4;i++)
	for(int t = 0;t < 3;t++)                                        
	    subver[i][t] = Xrotate[t][0] * rotateVertex[i][0] + Xrotate[t][1] * rotateVertex[i][1] + Xrotate[t][2] * rotateVertex[i][2];
    
    for(int i = 0;i < 4;i++)
	for(int t = 0;t < 3;t++)
	    rotateVertex[i][t] = Zrotate[t][0] * subver[i][0] + Zrotate[t][1] * subver[i][1] + Zrotate[t][2] * subver[i][2]; 

    double weightsum = 0;
    
    rotateVertex[0][0] += move;

    //t0,tを求める
    double t_0[3] = {0,0,0},t[3] = {0,0,0};      

    for(int i = 0;i < 4;i++) //点の数のループ
      {
	  weightsum += weight[i];
	  
	  for(int u = 0;u < 3;u++) //座標の数のループ　
	      {
		  t_0[u] += weight[i] * vertex[i][u];
		  t[u] += weight[i] * rotateVertex[i][u];
	      }
      }
    
    for(int i = 0;i < 3;i++)
	{
	    t_0[i] /= weightsum;
	    t[i] /= weightsum;
	}
  
    double p[3]; //変形後
    double q[3]; //変形前
    double ptq[3][3]; //p*q^t
    double qtq[3][3]; //q*q^t
    double sum_mptq[3][3]; //sum of m*ptq
    double sum_mqtq[3][3]; //sum of m*qtq
    
    //sumの初期化
    for(int i = 0;i < 3;i++)
      {
	  for(int t = 0;t < 3;t++)
	      {
		  sum_mqtq[i][t] = 0;
		  sum_mptq [i][t] = 0;
	      }
      }
  //sumの計算
  for(int i = 0;i < 4;i++)
      {
	  for(int u = 0;u < 3;u++)
	      {
		  q[u] = vertex[i][u] - t_0[u];
		  p[u] = rotateVertex[i][u] - t[u];
	      }
	  for(int u = 0;u < 3;u++)
	      {
		  for(int s = 0;s < 3;s++)
		      {
			  sum_mptq[u][s] += weight[i] * (p[u] * q[s]);
			  sum_mqtq[u][s] += weight[i] * (q[u] * q[s]);
		      }
	      }
      }

  double tsum_mptq[3][3]; //(sum_mptq)^t
  double inv_sum_mqtq[3][3]; //inverse of sum_mqtq
  double A[3][3]; 
  double S[3][3];

  for(int i = 0;i < 3;i++)
      {
	  for(int t = 0;t < 3;t++)
	      {
		  tsum_mptq[t][i] = sum_mptq[i][t];
	      }
      }

  Inverse(sum_mqtq,inv_sum_mqtq);

  Matrix *eig;
  Matrix *test;
  Vector *eigenvec;

  eig = CreateMatrix(3,3);
  test = CreateMatrix(3,3);

  //calc A and S
  for(int i = 0;i < 3;i++)
      {
	  for(int t = 0;t < 3;t++)
	      {
		  A[i][t] = sum_mptq[i][0] * inv_sum_mqtq[0][t] + sum_mptq[i][1] * inv_sum_mqtq[1][t] + sum_mptq[i][2] * inv_sum_mqtq[2][t];

		  test->a[i * test->width + t] = tsum_mptq[i][0] * sum_mptq[0][t] + tsum_mptq[i][1] * sum_mptq[1][t] + tsum_mptq[i][2] * sum_mptq[2][t];
	      }
      }
  
  eigenvec = Jacobi(test,eig);

  double P[3][3];
  double inv_P[3][3];
  double EIG[3][3];

  for(int i=0; i<eig->height; i++){
      for(int j=0; j<eig->width; j++){
	  P[i][j] = eig->a[i * eig->width + j];
      }
  }

  for(int i = 0; i < eigenvec->size;i++)
      {
          for(int t = 0;t < 3;t++)
              {
                  if(i == t) EIG[i][i] = sqrt(eigenvec->v[i]);
                  else EIG[i][t] = 0;
              }
      }

  //p*EIG*invp                                                                      
  double tmp[3][3];

  for(int i = 0;i < 3;i++)
      {
          for(int t = 0;t < 3;t++)
              {
		  tmp[i][t] = P[i][0] * EIG[0][t] + P[i][1] * EIG[1][t] + P[i][2] * EIG[2][t];
              }
      }

  Inverse(P,inv_P);

  for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
      }
  }

  for(int i = 0;i < 3;i++)                                                                   for(int t = 0;t < 3;t++)
	  S[i][t] = tmp[i][0] * inv_P[0][t] + tmp[i][1] * inv_P[1][t] + tmp[i][2] * inv_P[2][t];

  double inv_S[3][3];
  
  Inverse(S,inv_S);
  
  double R[3][3];

  for(int i = 0;i < 3;i++)
      for(int t = 0;t < 3;t++)
	  R[i][t] = tsum_mptq[0][i] * inv_S[0][t] + tsum_mptq[1][i] * inv_S[1][t] + tsum_mptq[2][i] * inv_S[2][t];
  
  for(int i = 0;i < 4;i++)
      {
	  for(int u = 0;u < 3;u++)
	      {
		  q[u] = vertex[i][u] - t_0[u];
	      }
	  for(int u = 0;u < 3;u++)                                                
              {   
		  g[i][u] = (R[u][0] * q[0] + R[u][1] * q[1] + R[u][2] * q[2]) + t[u];
	      }		  
      }

  //移動後の重心を求める
  double t_1[3] = {0,0,0};
                                                                                                          
  for(int i = 0;i < 4;i++) //点の数のループ                                                                                               
      for(int u = 0;u < 3;u++) //座標の数のループ　                                                               
	  t_1[u] += weight[i] * g[i][u];

  for(int i = 0;i < 3;i++) t_1[i] /= weightsum;  
  
  FILE *fp;

  if((fp = fopen("check-vertex.dat","w")) == NULL)
      {
	  printf("chck-vertex.dat can't open.\n");
	  exit(1);
      }
  for(int i = 0;i < 4;i++)
  fprintf(fp,"(%lf,%lf,%lf)\n(%lf,%lf,%lf)\n\n",rotateVertex[i][0],rotateVertex[i][1],rotateVertex[i][2],g[i][0],g[i][1],g[i][2]);

  fprintf(fp,"(%lf,%lf,%lf)\n(%lf,%lf,%lf)\n\n",t[0],t[1],t[2],t_1[0],t_1[1],t_1[2]);

  FreeVector(eigenvec);
  FreeMatrix(test);
  FreeMatrix(eig);
  fclose(fp);

  //線形有限粗放の合成マトリクスを計算
}

void idle(void)
{                                                  
   glutPostRedisplay();
 }

//軸を描画
void draw_axis()
{
    //x軸を描画
    glColor3d(1.0,0.0,0.0); //赤色
    glBegin(GL_LINES);
    glVertex3d(-100,0,0); //始点
    glVertex3d(100,0,0); //終点
    glEnd();
    //y軸を描画
    glColor3d(0.0,0.0,1.0); //青色
    glBegin(GL_LINES);
    glVertex3d(0,100,0); //始点
    glVertex3d(0,-100,0); //終点
    glEnd();
    //z軸を描画
    glColor3d(0.0,1.0,0.0); //緑色
    glBegin(GL_LINES);
    glVertex3d(0,0,100); //始点
    glVertex3d(0,0,-100); //終点
    glEnd();
}

//描画処理関数
void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT); //ウィンドウを背景色で塗りつぶし

  draw_axis();
  Shape_Matching();

  //原図形
  glPushMatrix();
  glColor3d(1,0,0);
  glBegin(GL_LINES);
  for(int i = 0;i < 6;i++)
      {
	  glVertex3dv(vertex[edge[i][0]]);
	  glVertex3dv(vertex[edge[i][1]]);
      }
  glEnd();
  glPopMatrix();
 
  glPushMatrix();
  glColor3d(1,0,0);
  glTranslated(vertex[0][0],vertex[0][1],vertex[0][2]);
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(0,1,0);
  glTranslated(vertex[1][0],vertex[1][1],vertex[1][2]);                           
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(0,0,1);
  glTranslated(vertex[2][0],vertex[2][1],vertex[2][2]);                           
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(1,0,1);
  glTranslated(vertex[3][0],vertex[3][1],vertex[3][2]);
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();         

   //回転図形  
  glPushMatrix();                                                        
  glColor3d(0,0,1);
  glBegin(GL_LINES);
  for(int i = 0;i < 6;i++)
      {
            glVertex3dv(rotateVertex[edge[i][0]]);
	    glVertex3dv(rotateVertex[edge[i][1]]);
      }
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(1,0,0);
  glTranslated(rotateVertex[0][0],rotateVertex[0][1],rotateVertex[0][2]);
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(0,1,0);
  glTranslated(rotateVertex[1][0],rotateVertex[1][1],rotateVertex[1][2]);
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(0,0,1);
  glTranslated(rotateVertex[2][0],rotateVertex[2][1],rotateVertex[2][2]);
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(1,0,1);
  glTranslated(rotateVertex[3][0],rotateVertex[3][1],rotateVertex[3][2]);
  glutSolidSphere(0.05,10,10);                                                    
  glEnd();
  glPopMatrix();

  //修正後
  glPushMatrix();                                                         
  glColor3d(0,1,0);
  glBegin(GL_LINES);
  for(int i = 0;i < 6;i++)
      {
	  glVertex3dv(g[edge[i][0]]);
	  glVertex3dv(g[edge[i][1]]);
      }
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(1,0,0);
  glTranslated(g[0][0],g[0][1],g[0][2]);
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(0,1,0);
  glTranslated(g[1][0],g[1][1],g[1][2]);
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(0,0,1);
  glTranslated(g[2][0],g[2][1],g[2][2]);
  glutSolidSphere(0.05,10,10);
  glEnd();
  glPopMatrix();

  glPushMatrix();
  glColor3d(1,0,1);
  glTranslated(g[3][0],g[3][1],g[3][2]);
  glutSolidSphere(0.05,10,10);                                                    
  glEnd();
  glPopMatrix();

  glutSwapBuffers();
}

//ウィンドウの拡大・縮小処理関数
void resize(int w, int h)
{
  /* ウィンドウ全体をビューポートにする */
  glViewport(0, 0, w, h);

  /* 透視変換行列の設定 */
  glMatrixMode(GL_PROJECTION);

  /* 変換行列の初期化 */
  glLoadIdentity();

  gluPerspective(30.0, (double)w / (double)h, 1.0, 100.0);

  glMatrixMode(GL_MODELVIEW);

  gluLookAt(3,3,7,0,0,0,0,1,0);
 }

void keyboard(unsigned char key, int x,int y)
{
  switch(key)
    {
    case 'x': r_x += 1.0; glutPostRedisplay(); break;
    case 'X': r_x -= 1.0; glutPostRedisplay(); break;
    case 'y': r_y += 1.0; glutPostRedisplay(); break;                                  
    case 'Y': r_y -= 1.0; glutPostRedisplay(); break;
    case 'z': r_z += 1.0; glutPostRedisplay(); break;                                   
    case 'Z': r_z -= 1.0; glutPostRedisplay(); break;                   
    case 'm': move += 1.0; glutPostRedisplay(); break;
    case 'M': move -= 1.0; glutPostRedisplay(); break;
    case 'q':
    case 'Q': exit(0);
    default: break;
    }
}

//ウィンドウの初期化関数
void init(void)
{
  /* ウィンドウの背景色の設定 */
  glClearColor(1.0, 1.0, 1.0, 1.0);
}
  
//掃き出し法                                                                                            
void Inverse(double mat[][3],double inv[][3])
{
    double buf;

    //invの初期化                                                                                       
    for(int i = 0;i < 3;i++)
	{
	    for(int t = 0;t < 3;t++)
		{
		    inv[i][t] = (i == t)?1.0:0.0;
		}
	}
    
        for(int i = 0;i < 3;i++)
	{
	    buf = 1 / mat[i][i];
	    for(int j = 0;j < 3;j++)
		{
		    mat[i][j] *= buf;
		    inv[i][j] *= buf;
		}
	    for(int t = 0;t < 3;t++)
		{
		    if(i != t)
			{
			    buf = mat[t][i];
			    for(int j = 0;j < 3;j++)
				{
				    mat[t][j] -= (buf * mat[i][j]);
				    inv[t][j] -= (buf * inv[i][j]);
				}
			}
		}
	}
    
}    


int main(int argc, char *argv[])
{
  cout << "Start OpenGL Program" << endl;
  glutInitWindowPosition(100, 100);  // ウィンドウの位置を指定
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);      // ウィンドウの幅と高さを指定
  glutInit(&argc, argv);
  glutIdleFunc(idle);
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE);
  glutCreateWindow(argv[0]);
  glutDisplayFunc(display);
  glutReshapeFunc(resize);
  glutKeyboardFunc(keyboard);
  init();
  glutMainLoop();
  return 0;
}
