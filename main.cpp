
/*
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include <windows.h>
#include <glut.h>
#include<queue>

*/
#include <iostream>
#include <fstream>
#include "bitmap_image.hpp"
#include "1605114_classes.h"
#define pi (2*acos(0.0))

using namespace std ;

vector <Object*> objects;
vector <Light> lights;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double ux,uy,uz,rx,ry,rz,lx,ly,lz,posx,posy,posz ;
double q_w_angle ;
double e_r_angle;
double a_s_angle ;
double d_f_angle ;
double board_distance ;

int  level_of_recursion ;
int pixels ;
int number_of_objects , number_of_lights ;

double viewAngle = 80 ;
double aspectRatio = 1 ;
double nearDistance = 1 ;
double farDistance = 1000;

int windowHeight = 500 ;
int windowWidth = 500 ;

//extern vector <Object*> objects;
//extern vector <Light> lights;

point pos2, u2,l2,r2 ;
point q[2000];
int queue_size;



void Capture()
{
    int pixels = 768 ;
    cout << "inside capture " << endl;
    Vector3D eye(posx,posy,posz);

    Color** image_color ;
    image_color = new Color*[pixels];

    for(int i = 0; i < pixels; ++i) {
        image_color[i] = new Color[pixels];
        for (int j = 0; j < pixels ; ++j) {
            Color color(0,0,0);
            image_color[i][j] = color ;
        }
    }

    double planeDistance = (windowHeight/2.0) / tan(viewAngle*(3.1416/180)*0.5);
    Vector3D topleft ;
    topleft.x = posx + lx*planeDistance - rx*windowWidth/2 +
                                    ux*windowHeight/2 ;
    topleft.y = posy + ly*planeDistance - ry*windowWidth/2 +
                                    uy*windowHeight/2 ;
    topleft.z = posz + lz*planeDistance - rz*windowWidth/2 +
                                    uz*windowHeight/2 ;


    double du = (windowWidth*1.0)/(pixels*1.0) ;
    double dv = (windowHeight*1.0)/(pixels*1.0)  ;

    topleft.x = topleft.x + rx*(0.5*du) - ux*(0.5*dv) ;
    topleft.y = topleft.y + ry*(0.5*du) - uy*(0.5*dv) ;
    topleft.z = topleft.z + rz*(0.5*du) - uz*(0.5*dv) ;

    cout <<"rx ----->  "<<rx << "  ux----> " << ux << " "<< endl;
    cout << "du --> "<< du << "  dv :  " << dv << endl;



    for(int i = 0 ; i < pixels ; i++){
        for(int j = 0 ; j < pixels ; j++){
            int nearest;
            double t , tMin = -1;

            Vector3D currentPixel ;
            currentPixel.x = topleft.x + rx*(du)*i - ux*(dv)*j ;
            currentPixel.y = topleft.y + ry*(du)*i - uy*(dv)*j ;
            currentPixel.z = topleft.z + rz*(du)*i - uz*(dv)*j ;

            Ray r(eye,currentPixel);
            double *color = new double[3];
            color[0] = 0; color[1] = 0 ; color[2] = 0 ;
            Object *o ;
            for(int k = 0 ; k < objects.size() ; k++){

                t = objects[k]->intersect(r,color,0);

                if(tMin < 0 && t > 0){

                    tMin = t ;
                    o = objects[k];

                }else if(tMin > t && t > 0){
                    tMin = t ;
                    o = objects[k];
                }

            }
            if(tMin > 0){
                tMin = o->intersect(r,color,1);
                Color c(color[0],color[1],color[2]);
                image_color[i][j] = c;
            }

            //baki.........


        }

    }



       bitmap_image image(pixels,pixels);

        for (int i = 0; i < pixels ; ++i) {
            for (int j = 0; j < pixels ; ++j) {
                image.set_pixel(i,j,image_color[i][j].color[0]*255,
                                image_color[i][j].color[1]*255,
                                image_color[i][j].color[2]*255);

            }
        }


    image.save_image("out.bmp");

    for( int i = 0 ; i < pixels ; i++ ){
        delete[] image_color[i];
    }

    delete[] image_color;

    cout << "Capture---------------------------------------------"<< endl;

}

void drawAxes()
{
	if(drawaxes==1)
	{
		//glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{

			glColor3f(1.0, 0, 0);             //red --X //green --> y blue --> z
			glVertex3f( 1000,0,0);
			glVertex3f(0,0,0);

			glColor3f(1, 1, 1.0);
			glVertex3f(0,0,0);
			glVertex3f(-1000,0,0);

			glColor3f(0, 1, 0);
			glVertex3f(0,1000,0);
			glVertex3f(0,0,0);

			glColor3f(1, 1, 1.0);
			glVertex3f(0,0,0);
			glVertex3f(0, -1000,0);

			glColor3f(0, 0, 1.0);
			glVertex3f(0,0, 1000);
			glVertex3f(0,0, 0);

			glColor3f(1,1,1.0);
			glVertex3f(0,0, 0);
			glVertex3f(0,0,-1000);

		}glEnd();
	}
}

void drawSS()
{

    //cout << "drawss : \n" << endl;
    for(int i = 0 ; i < objects.size() ; i++)
    {
        objects[i]->draw();
    }
    for(int i = 0 ; i < lights.size() ; i++)
    {
        lights[i].draw();
    }


    //objects[0]->draw();
}

void keyboardListener(unsigned char key, int x,int y){
	double theta = 0.05 ;
	double sin_theta = sin(theta);
	double cos_theta = cos(theta);
	//cout << sin_theta << " " << cos_theta << endl;
	switch(key){
        case '0':
        {
            Capture();
            break;
        }
		case '1':
		    {
		        double prevlx = lx;
                double prevly = ly;
                double prevlz = lz;
                //drawgrid=1-drawgrid;
                lx = -rx* sin_theta + lx*cos_theta;
                ly = -ry* sin_theta + ly*cos_theta;
                lz = -rz* sin_theta + lz*cos_theta;

                rx = prevlx*sin_theta + rx*cos_theta;
                ry = prevly*sin_theta + ry*cos_theta;
                rz = prevlz*sin_theta + rz*cos_theta;

                //cout << lx << " " << ly << " "<< lz << endl;
                //cout << rx << " " << ry << "  " << rz << endl;
                break;
		    }
        case '2':
		    {
		        double prevlx = lx;
                double prevly = ly;
                double prevlz = lz;
                //drawgrid=1-drawgrid;
                lx = rx* sin_theta + lx*cos_theta;
                ly = ry* sin_theta + ly*cos_theta;
                lz = rz* sin_theta + lz*cos_theta;

                rx = -prevlx*sin_theta + rx*cos_theta;
                ry = -prevly*sin_theta + ry*cos_theta;
                rz = -prevlz*sin_theta + rz*cos_theta;

                //cout << lx << " " << ly << " "<< lz << endl;
                //cout << rx << " " << ry << "  " << rz << endl;
                break;
		    }
        case '3':
		    {
		        double prevlx = lx;
                double prevly = ly;
                double prevlz = lz;
                //drawgrid=1-drawgrid;
                lx = ux* sin_theta + lx*cos_theta;
                ly = uy* sin_theta + ly*cos_theta;
                lz = uz* sin_theta + lz*cos_theta;

                ux = -prevlx*sin_theta + ux*cos_theta;
                uy = -prevly*sin_theta + uy*cos_theta;
                uz = -prevlz*sin_theta + uz*cos_theta;

                break;
		    }
        case '4':
		    {
		        double prevlx = lx;
                double prevly = ly;
                double prevlz = lz;
                //drawgrid=1-drawgrid;
                lx = -ux* sin_theta + lx*cos_theta;
                ly = -uy* sin_theta + ly*cos_theta;
                lz = -uz* sin_theta + lz*cos_theta;

                ux = prevlx*sin_theta + ux*cos_theta;
                uy = prevly*sin_theta + uy*cos_theta;
                uz = prevlz*sin_theta + uz*cos_theta;

                break;
		    }
        case '5':
		    {
		        double prevUx = ux;
                double prevUy = uy;
                double prevUz = uz;
                //drawgrid=1-drawgrid;

                ux = -rx*sin_theta + ux*cos_theta;
                uy = -ry*sin_theta + uy*cos_theta;
                uz = -rz*sin_theta + uz*cos_theta;

                rx = prevUx* sin_theta + rx*cos_theta;
                ry = prevUy* sin_theta + ry*cos_theta;
                rz = prevUz* sin_theta + rz*cos_theta;

                break;
		    }
        case '6':
		    {

		        double prevUx = ux;
                double prevUy = uy;
                double prevUz = uz;
                //drawgrid=1-drawgrid;

                ux = rx*sin_theta + ux*cos_theta;
                uy = ry*sin_theta + uy*cos_theta;
                uz = rz*sin_theta + uz*cos_theta;

                rx = -prevUx* sin_theta + rx*cos_theta;
                ry = -prevUy* sin_theta + ry*cos_theta;
                rz = -prevUz* sin_theta + rz*cos_theta;


                break;
		    }

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			posx -=10*lx;
			posy -=10*ly;
			posz -=10*lz;
			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;
			posx +=10*lx;
			posy +=10*ly;
			posz +=10*lz;
			break;

		case GLUT_KEY_RIGHT:
			posx += 10*rx;
			posy += 10*ry;
			posz += 10*rz;
			//posy +=2*posy;
			//posz +=2*posz;
			break;
		case GLUT_KEY_LEFT:
			posx -= 10*rx;
			posy -= 10*ry;
			posz -= 10*rz;
			break;

		case GLUT_KEY_PAGE_UP:
		    posx += 10*ux;
			posy += 10*uy;
			posz += 10*uz;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    posx -= 10*ux;
			posy -= 10*uy;
			posz -= 10*uz;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP


			}
			break;

		case GLUT_RIGHT_BUTTON:
		    if(state == GLUT_DOWN){
                drawaxes=1-drawaxes;
		    }
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
	//cout << "In display : " << posx+lx << "  " << posy+ly << "  "<< posz+lz <<endl;
    gluLookAt(posx,posy,posz,	posx+lx,posy+ly,posz+lz,	ux,uy,uz);

    //posx,posy,posz __> camera

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
    drawSS();

	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
	ux = 0,uy=0,uz=1;
	rx= (-1/sqrt(2)),ry=(1/sqrt(2)),rz=0;
	lx=(-1/sqrt(2)),ly=(-1/sqrt(2)),lz=0;
	posx=100,posy=100,posz=0 ;
	q_w_angle = 0;
    e_r_angle = 0;
    a_s_angle = 0;
    d_f_angle = 0;
    queue_size = 0;
    board_distance = 300;

    //point pos2, u2,l2,r2 ;
    /*
    pos2.x = 0; pos2.y = 40; pos2.z = 0;
    l2.x = 0; l2.y = 340; l2.z = 0;
    u2.x = 0; u2.y = 40; u2.z = 300;
    r2.x = 300; r2.y = 40; r2.z = 0;
    */
	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters

	gluPerspective(viewAngle,	aspectRatio ,nearDistance,	farDistance);  // view angle --> 80
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}




void loadData()
{

    fstream scenefile;
    scenefile.open("scene.txt");
    if (scenefile.is_open()) {
        scenefile >> level_of_recursion >> pixels >> number_of_objects ;
        cout << level_of_recursion << " " << pixels << endl ;

        int it_count = 0 ;
        string input_command ;
        while (it_count < number_of_objects){
            it_count++ ;
            scenefile >> input_command ;
            if (input_command == "sphere"){

                double x , y , z , radius, r,g,b ,ambient, diffuse,
                        specular, recursive_reflection_coefficient, shininess  ;
                scenefile >> x >> y >> z >> radius ;
                scenefile >> r >> g >> b ;
                scenefile >> ambient >> diffuse >> specular >>
                        recursive_reflection_coefficient >> shininess ;
                Vector3D center(x,y,z);
                Object *temp = new Sphere(center,radius);
                temp->setColor(r,g,b);
                temp->setShine(shininess);
                temp->setCoEfficients(ambient,diffuse,specular,recursive_reflection_coefficient);
                objects.push_back(temp);

            }else if(input_command == "triangle"){

                double x1 , y1 , z1 , x2 , y2 , z2 ,x3 , y3 , z3 , r , g , b ,ambient, diffuse,
                        specular, recursive_reflection_coefficient, shininess  ;
                scenefile >> x1 >> y1 >> z1 >>  x2 >> y2 >> z2 >>  x3 >> y3 >> z3  ;
                scenefile >> r >> g >> b ;
                scenefile >> ambient >> diffuse >> specular >>
                        recursive_reflection_coefficient >> shininess ;
                Vector3D p1(x1,y1,z1);
                Vector3D p2(x2,y2,z2);
                Vector3D p3(x3,y3,z3);
                Object *temp = new Triangle(p1,p2,p3);
                temp->setColor(r,g,b);
                temp->setShine(shininess);
                temp->setCoEfficients(ambient,diffuse,specular,recursive_reflection_coefficient);
                objects.push_back(temp);


            }else if(input_command == "general"){

                double A ,B ,C ,D ,E ,F ,G ,H ,I ,J ,length , width, height, x , y , z , r , g , b ,ambient, diffuse,
                        specular, recursive_reflection_coefficient, shininess  ;

                scenefile >> A >> B >> C >> D >>E >> F >> G >> H >> I >> J  ;
                scenefile >> x >> y >> z  ;
                scenefile >> length >> width >> height  ;
                scenefile >> r >> g >> b ;
                scenefile >> ambient >> diffuse >> specular >>
                        recursive_reflection_coefficient >> shininess ;


                Vector3D ref_point(x,y,z);

                Object *temp = new General(A,B,C,D,E,F,G,H,I,J,ref_point,length,width,height);
                temp->setColor(r,g,b);
                temp->setShine(shininess);
                temp->setCoEfficients(ambient,diffuse,specular,recursive_reflection_coefficient);
                objects.push_back(temp);

            }

        }
        scenefile >>  number_of_lights ;
        it_count = 0 ;
        while(it_count < number_of_lights)
        {
            it_count++ ;
            double x , y , z , r ,g ,b;
            scenefile >> x >> y >> z ;
            scenefile >> r >> g >> b ;
            Vector3D position(x,y,z) ;
            Light light(position, r,g,b);
            lights.push_back(light);
        }

        scenefile.close();
    }


    Object *temp = new Floor(1000,20);
    temp->setColor(1,1,1);
    temp->setCoEfficients(0.2,0.2,0.3,0.4);
    temp->setShine(10);
    objects.push_back(temp);


}


int main(int argc, char **argv){
    //test();
    loadData();

	glutInit(&argc,argv);
	glutInitWindowSize(windowWidth,windowHeight);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

    //Capture();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	for (Object* obj : objects)
        delete obj;
    objects.clear();

	return 0;
}
