#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include <windows.h>
#include <glut.h>
#include<queue>
#define pi (2*acos(0.0))

using namespace std ;
class Object ;
class Light ;
extern vector <Object*> objects;
extern vector <Light> lights;
extern int  level_of_recursion ;

double sml = 0.0000001 ;
void test()
{
    cout << "from the header class "<< endl;
}



class Color
{
    public:
        double color[3];
        Color(){}
        Color(double r, double g, double b)
        {
            color[0] = r ;
            color[1] = g ;
            color[2] = b ;
        }

        void print(){
            cout << "color : " << "r : " << color[0]
            << "g : " << color[1]<< " b : "<< color[2] << endl;
        }
};

class Vector3D{
    public:
        double x , y , z ;
        Vector3D(double x1 , double y1 , double z1)
        {
            x = x1 ; y = y1 ; z = z1 ;
        }
        Vector3D()
        {
        }
        void normalize()
        {
            double d = sqrt(x*x + y*y + z*z) ;
            x = x/d ; y = y/d ; z = z/d ;
        }
        void print()
        {
            cout << x << " " << y << " " << z << endl;
        }
        double magnitude()
        {
            return sqrt(x*x+y*y+z*z);
        }
};

struct point
{
	double x,y,z;
};

double dot_product(Vector3D a, Vector3D b)
{
    double product = 0;
    product += a.x * b.x;
    product += a.y * b.y;
    product += a.z * b.z;
    return product;
}

Vector3D cross_product(Vector3D a, Vector3D b)
{
    Vector3D result ;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

Vector3D subtract_vector(Vector3D a, Vector3D b)
{
    Vector3D result( a.x-b.x , a.y-b.y , a.z-b.z );
    return result ;
}

Vector3D scaleVector(Vector3D a, double val)
{
    Vector3D result ;
    result.x = a.x * val ;
    result.y = a.y * val;
    result.z = a.z * val;
    return result;
}

class Matrix
{
    public :
        double matrix[3][3];

        double determinant( double matrix[3][3], int n) {
           double det = 0;
           double submatrix[3][3];
           if (n == 2)
           return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
           else {
              for (int x = 0; x < n; x++) {
                 int subi = 0;
                 for (int i = 1; i < n; i++) {
                    int subj = 0;
                    for (int j = 0; j < n; j++) {
                       if (j == x)
                       continue;
                       submatrix[subi][subj] = matrix[i][j];
                       subj++;
                    }
                    subi++;
                 }
                 det = det + (pow(-1, x) * matrix[0][x] * determinant( submatrix, n - 1 ));
              }
           }
           return det;
        }

};



class Ray{
    public:
        Vector3D start;
        Vector3D dir; // normalize for easier calculations
        //write appropriate constructor
        Ray(Vector3D start1, Vector3D end1 )
        {
            start = start1;
            Vector3D direction = subtract_vector(end1,start1);
            direction.normalize();
            dir = direction ;
        }
        Ray(){}
};

class Light{
    public:
        Vector3D light_pos;
        double color[3];

        Light(Vector3D light_pos_ , double red , double green , double blue )
        {
            //cout << "Light Constructor " << endl;
            light_pos = light_pos_ ;
            color[0] = red ;
            color[1] = green ;
            color[2] = blue ;
        }
        void draw(){
            //cout << "In Light draw : " << endl;
            glColor3f(color[0],color[1],color[2]);
            glPushMatrix();
            glTranslatef(light_pos.x,light_pos.y,light_pos.z);
            glutSolidSphere(2,100,100);
            glPopMatrix();
        }
};

class Object{
    public:
        Vector3D reference_point ;// should have x, y, z
        double height, width, length ;
        double color[3] ;
        double coEfficients[4] ;// reflection coefficients
        int shine ; // exponent term of specular component

        virtual ~Object(){
            cout << "Object constructor" << endl;
        }
        virtual void draw(){
            cout << "In object draw : " << endl;
        }
        void setColor(double r, double g , double b){
            color[0] = r ;
            color[1] = g ;
            color[2] = b ;
        }
        void setShine(int s){
            shine = s ;
        }
        void setCoEfficients(double ambient,double diffuse, double specular,double recursive_reflection_coefficient){
            coEfficients[0] = ambient ;
            coEfficients[1] = diffuse ;
            coEfficients[2] = specular ;
            coEfficients[3] = recursive_reflection_coefficient ;

        }
        virtual double intersect(Ray r, double *col, int level){
            return -1.0;
        }
};



class Sphere : public Object{
    public:
        Sphere(Vector3D center, double radius){
            //cout << "Sphere Constructor "<< endl;
            reference_point = center ;
            length = radius ;
        }
        void draw(){
            //cout << "In sphere draw : " << endl;
            glColor3f(color[0],color[1],color[2]);
            glPushMatrix();
            glTranslatef(reference_point.x,reference_point.y,reference_point.z);
            glutSolidSphere(length,100,100);
            glPopMatrix();
        }
        double intersect(Ray r, double *col, int level){
            //cout << "sphere inresect "<< endl;

            Vector3D Ro = subtract_vector(r.start,reference_point ); //check
            double a = 1 ;
            double b = 2* dot_product(r.dir, Ro);
            double c = dot_product(Ro,Ro) - length*length ;


            double d_square = b*b - 4*a*c ;
            if(d_square < 0){
                return -1 ;
            }
            double d = sqrt(d_square) ;
            double t1 = (-b+d)/(2*a);
            double t2 = (-b-d)/(2*a);



            double t ;
            if(t1< 0 && t2 < 0){
                t =  -1 ;
            }else if(t1<0 && t2>0){
                t = t2 ;
            }else if(t1>0 && t2<0){
                t = t1 ;
            }else if(t1>0 && t2>0 && t1>t2 ){
                t = t2 ;
            }else if(t1>0 && t2>0 && t1<t2 ){
                t = t1 ;
            }

            if(level == 0){
                return t ;
            }



            Vector3D intersectionPoint ;
            intersectionPoint.x = r.start.x + t*r.dir.x ;
            intersectionPoint.y = r.start.y + t*r.dir.y ;
            intersectionPoint.z = r.start.z + t*r.dir.z ;

            double intersectionPointColor[3];
            intersectionPointColor[0] = color[0];
            intersectionPointColor[1] = color[1];
            intersectionPointColor[2] = color[2];

            col[0] = intersectionPointColor[0]*coEfficients[0] ;
            col[1] = intersectionPointColor[1]*coEfficients[0] ;
            col[2] = intersectionPointColor[2]*coEfficients[0] ;


            Vector3D normal = subtract_vector(intersectionPoint,reference_point);
            normal.normalize();

            for(int i = 0 ; i < lights.size() ; i++){
                Ray ray1(intersectionPoint,lights[i].light_pos);
                ray1.start = Vector3D(intersectionPoint.x+ sml*ray1.dir.x ,
                                      intersectionPoint.y+ sml*ray1.dir.y ,
                                      intersectionPoint.z+ sml*ray1.dir.z);


                int obscured = 0 ;
                double tMin = -1 ;
                for(int k = 0 ; k < objects.size() ; k++){
                    t = objects[k]->intersect(ray1,color,0);
                    if(tMin < 0 && t > 0){
                        tMin = t ;
                    }else if(tMin > t && t > 0){
                        tMin = t ;
                    }
                }
                Vector3D distance_vec = subtract_vector(lights[i].light_pos, intersectionPoint);
                double distance = distance_vec.magnitude();

                if(tMin > 0 && distance > tMin){
                    obscured = 1 ;
                }

                if(obscured == 0 )
                {
                    //cout << "Here" << endl;
                    double lambert = dot_product(normal,ray1.dir) > 0 ?
                                      dot_product(normal,ray1.dir) : 0;
                    //find reflected ray, rayr for rayl
                    //R = 2 (L.N)N – L
                    Ray ray_r ;
                    ray_r.start = intersectionPoint ;
                    ray_r.dir = subtract_vector(scaleVector(normal,2*lambert),ray1.dir);
                    //ray_r.dir = subtract_vector(ray1.dir,scaleVector(normal,2*lambert));
                    ray_r.dir.normalize() ;

                    // r k negative kore
                    Vector3D negative_r = scaleVector(r.dir,-1);
                    negative_r.normalize();
                    double phongValue =  dot_product(negative_r,ray_r.dir) > 0 ?
                                dot_product(negative_r,ray_r.dir) : 0;

                    //cout << "lambert : "<< lambert << "  phong : " << phongValue << endl;
                    col[0] += lights[i].color[0]*coEfficients[1]*lambert*intersectionPointColor[0] ;
                    col[0] += lights[i].color[0]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[0] ;

                    col[1] += lights[i].color[1]*coEfficients[1]*lambert*intersectionPointColor[1] ;
                    col[1] += lights[i].color[1]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[1] ;

                    col[2] += lights[i].color[2]*coEfficients[1]*lambert*intersectionPointColor[2] ;
                    col[2] += lights[i].color[2]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[2] ;

                     if(col[0] < 0){
                            col[0] = 0 ;
                        }
                        if(col[1] < 0){
                            col[1] = 0 ;
                        }
                        if(col[2] < 0){
                            col[2] = 0 ;
                        }
                        if(col[0] > 1){
                            col[0] = 1 ;
                        }
                        if(col[1] > 1){
                            col[1] = 1 ;
                        }
                        if(col[2] > 1){
                            col[2] = 1 ;
                        }



                }
            }


            if(level >= level_of_recursion ){
                //cout << "level 4 : "<< col[0] << " "<< col[1] << " " << col[2] << endl;
                return t ;
            }

            Ray reflected ;
            //R = 2 (L.N)N – L
            double scale = 2* dot_product(r.dir,normal);
            //reflected.dir = subtract_vector(scaleVector(normal,scale),r.dir);
            reflected.dir = subtract_vector(r.dir , scaleVector(normal,scale));

            reflected.dir.normalize();
            reflected.start = Vector3D(intersectionPoint.x + sml*reflected.dir.x ,
                              intersectionPoint.y+ sml*reflected.dir.y ,
                              intersectionPoint.z+ sml*reflected.dir.z);

            double t_ , tMin_ = -1 ;
            double *reflectedColor = new double[3];
            reflectedColor[0] = 0.0; reflectedColor[1] = 0.0 ; reflectedColor[2] = 0.0 ;
            Object *o ;
            for(int k = 0 ; k < objects.size() ; k++){
                t_ = objects[k]->intersect(reflected , reflectedColor,0);
                if(tMin_ < 0 && t_ > 0){
                    tMin_ = t_ ;
                    o = objects[k];
                }else if(tMin_ > t_ && t_ > 0){
                    tMin_ = t_ ;
                    o = objects[k];
                }
            }
            if(tMin_ > 0){
                o->intersect(reflected, reflectedColor ,level+1);
                //cout << "level : " << level << endl;

                //cout << reflectedColor[0] << " " << reflectedColor[1] << " " << reflectedColor[2] << endl;
                col[0] += reflectedColor[0] * coEfficients[3];
                col[1] += reflectedColor[1] * coEfficients[3];
                col[2] += reflectedColor[2] * coEfficients[3];
            }




            if(col[0] < 0){
                col[0] = 0 ;
            }
            if(col[1] < 0){
                col[1] = 0 ;
            }
            if(col[2] < 0){
                col[2] = 0 ;
            }
            if(col[0] > 1){
                col[0] = 1 ;
            }
            if(col[1] > 1){
                col[1] = 1 ;
            }
            if(col[2] > 1){
                col[2] = 1 ;
            }

            //return t ;
            //cout << "MMMMMMMMMMMMMMMM" << endl;
        }
};

class Triangle : public Object{

    public:
        Vector3D p1,p2,p3;

        Triangle(Vector3D p_1, Vector3D p_2, Vector3D p_3){
            //cout << "Triangle Constructor "<< endl;
            p1 = p_1;
            p2 = p_2;
            p3 = p_3;
        }
        void draw(){
            //cout << "In triangle draw : " << endl;
            glBegin(GL_TRIANGLES);{
                glColor3f(color[0],color[1],color[2]);
                glVertex3f( p1.x, p1.y,p1.z);
                glVertex3f( p2.x, p2.y,p2.z);
                glVertex3f( p3.x, p3.y,p3.z);
                //glVertex3f(-a, a,2);
            }glEnd();
        }

        double intersect(Ray r, double *col, int level){
            //cout << "Triangle inresect "<< endl;
            Matrix A , B ,G, T;
            A.matrix[0][0] = p1.x - p2.x ;  A.matrix[0][1] = p1.x - p3.x ; A.matrix[0][2] = r.dir.x ;
            A.matrix[1][0] = p1.y - p2.y ;  A.matrix[1][1] = p1.y - p3.y ; A.matrix[1][2] = r.dir.y ;
            A.matrix[2][0] = p1.z - p2.z ;  A.matrix[2][1] = p1.z - p3.z ; A.matrix[2][2] = r.dir.z ;

            B.matrix[0][0] = p1.x - r.start.x ;  B.matrix[0][1] = p1.x - p3.x ; B.matrix[0][2] = r.dir.x ;
            B.matrix[1][0] = p1.y - r.start.y ;  B.matrix[1][1] = p1.y - p3.y ; B.matrix[1][2] = r.dir.y ;
            B.matrix[2][0] = p1.z - r.start.z ;  B.matrix[2][1] = p1.z - p3.z ; B.matrix[2][2] = r.dir.z ;

            G.matrix[0][0] = p1.x - p2.x ;  G.matrix[0][1] = p1.x - r.start.x ; G.matrix[0][2] = r.dir.x ;
            G.matrix[1][0] = p1.y - p2.y ;  G.matrix[1][1] = p1.y - r.start.y ; G.matrix[1][2] = r.dir.y ;
            G.matrix[2][0] = p1.z - p2.z ;  G.matrix[2][1] = p1.z - r.start.z ; G.matrix[2][2] = r.dir.z ;

            T.matrix[0][0] = p1.x - p2.x ;  T.matrix[0][1] = p1.x - p3.x ; T.matrix[0][2] = p1.x - r.start.x ;
            T.matrix[1][0] = p1.y - p2.y ;  T.matrix[1][1] = p1.y - p3.y ; T.matrix[1][2] = p1.y - r.start.y ;
            T.matrix[2][0] = p1.z - p2.z ;  T.matrix[2][1] = p1.z - p3.z ; T.matrix[2][2] = p1.z - r.start.z ;

            double beta, gama , t ;
            beta = B.determinant(B.matrix,3)/ A.determinant(A.matrix,3) ;
            gama = G.determinant(G.matrix,3)/ A.determinant(A.matrix,3) ;
            t = T.determinant(T.matrix,3)/ A.determinant(A.matrix,3) ;

            //cout << "TRiagnle : " <<"t : "<<t <<"   Beta : " << beta << "  Gama : "<< gama <<  endl;
            if(beta+gama < 1 && beta > 0 && gama > 0 && t > 0 ){
                if(level == 0){
                    return t ;
                }


                //phong

                Vector3D intersectionPoint ;
                intersectionPoint.x = r.start.x + t*r.dir.x ;
                intersectionPoint.y = r.start.y + t*r.dir.y ;
                intersectionPoint.z = r.start.z + t*r.dir.z ;

                double intersectionPointColor[3];
                intersectionPointColor[0] = color[0];
                intersectionPointColor[1] = color[1];
                intersectionPointColor[2] = color[2];

                col[0] = intersectionPointColor[0]*coEfficients[0] ;
                col[1] = intersectionPointColor[1]*coEfficients[0] ;
                col[2] = intersectionPointColor[2]*coEfficients[0] ;

                Vector3D b_a = subtract_vector(p2,p1);
                Vector3D c_a = subtract_vector(p3,p1);
                Vector3D normal = cross_product(b_a,c_a);
                normal.normalize();
                normal.normalize();

                for(int i = 0 ; i < lights.size() ; i++){
                    Ray ray1(intersectionPoint,lights[i].light_pos);
                    ray1.start = Vector3D(intersectionPoint.x+ sml*ray1.dir.x ,
                                      intersectionPoint.y+ sml*ray1.dir.y ,
                                      intersectionPoint.z+ sml*ray1.dir.z);//Vector3D ray1(subtract_vector(intersectionPoint,lights[i].light_pos));


                    int obscured = 0 ;
                    double tMin = -1 ;
                    for(int k = 0 ; k < objects.size() ; k++){
                        t = objects[k]->intersect(ray1,color,0);
                        if(tMin < 0 && t > 0){
                            tMin = t ;
                        }else if(tMin > t && t > 0){
                            tMin = t ;
                        }
                    }
                    Vector3D distance_vec = subtract_vector(lights[i].light_pos, intersectionPoint);
                    double distance = distance_vec.magnitude();

                    if(tMin > 0 && distance > tMin){
                        obscured = 1 ;
                    }

                    //if(tMin > 0 && abs(t_this - tMin) < 0.00000001 )
                    if(obscured == 0 )
                    {
                        //cout << "Here" << endl;
                        double lambert = dot_product(normal,ray1.dir) > 0 ?
                                          dot_product(normal,ray1.dir) : 0;
                        //find reflected ray, rayr for rayl
                        //R = 2 (L.N)N – L
                        Ray ray_r ;
                        ray_r.start = intersectionPoint ;
                        //confused
                        ray_r.dir = subtract_vector(scaleVector(normal,2*lambert),ray1.dir);
                        //ray_r.dir = subtract_vector(ray1.dir,scaleVector(normal,2*lambert));
                        ray_r.dir.normalize() ;

                        // r k negative kore
                        Vector3D negative_r = scaleVector(r.dir,-1);
                        negative_r.normalize();
                        double phongValue =  dot_product(negative_r,ray_r.dir) > 0 ?
                                    dot_product(negative_r,ray_r.dir) : 0;

                        //cout << "lambert : "<< lambert << "  phong : " << phongValue << endl;
                        col[0] += lights[i].color[0]*coEfficients[1]*lambert*intersectionPointColor[0] ;
                        col[0] += lights[i].color[0]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[0] ;

                        col[1] += lights[i].color[1]*coEfficients[1]*lambert*intersectionPointColor[1] ;
                        col[1] += lights[i].color[1]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[1] ;

                        col[2] += lights[i].color[2]*coEfficients[1]*lambert*intersectionPointColor[2] ;
                        col[2] += lights[i].color[2]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[2] ;

                        if(col[0] < 0){
                            col[0] = 0 ;
                        }
                        if(col[1] < 0){
                            col[1] = 0 ;
                        }
                        if(col[2] < 0){
                            col[2] = 0 ;
                        }
                        if(col[0] > 1){
                            col[0] = 1 ;
                        }
                        if(col[1] > 1){
                            col[1] = 1 ;
                        }
                        if(col[2] > 1){
                            col[2] = 1 ;
                        }

                    }
            }

            if(level >= level_of_recursion ){
                //cout << "level 4 : "<< col[0] << " "<< col[1] << " " << col[2] << endl;
                return t ;
            }

            Ray reflected ;
            //R = 2 (L.N)N – L
            double scale = 2* dot_product(r.dir,normal);
            //reflected.dir = subtract_vector(scaleVector(normal,scale),r.dir);
            reflected.dir = subtract_vector(r.dir , scaleVector(normal,scale));

            reflected.dir.normalize();
            reflected.start = Vector3D(intersectionPoint.x + sml*reflected.dir.x ,
                              intersectionPoint.y+ sml*reflected.dir.y ,
                              intersectionPoint.z+ sml*reflected.dir.z);

            double t_ , tMin_ = -1 ;
            double *reflectedColor = new double[3];
            reflectedColor[0] = 0.0; reflectedColor[1] = 0.0 ; reflectedColor[2] = 0.0 ;
            Object *o ;
            for(int k = 0 ; k < objects.size() ; k++){
                t_ = objects[k]->intersect(reflected , reflectedColor,0);
                if(tMin_ < 0 && t_ > 0){
                    tMin_ = t_ ;
                    o = objects[k];
                }else if(tMin_ > t_ && t_ > 0){
                    tMin_ = t_ ;
                    o = objects[k];
                }
            }
            if(tMin_ > 0){
                o->intersect(reflected, reflectedColor ,level+1);
                col[0] += reflectedColor[0] * coEfficients[3];
                col[1] += reflectedColor[1] * coEfficients[3];
                col[2] += reflectedColor[2] * coEfficients[3];
            }




            if(col[0] < 0){
                col[0] = 0 ;
            }
            if(col[1] < 0){
                col[1] = 0 ;
            }
            if(col[2] < 0){
                col[2] = 0 ;
            }
            if(col[0] > 1){
                col[0] = 1 ;
            }
            if(col[1] > 1){
                col[1] = 1 ;
            }
            if(col[2] > 1){
                col[2] = 1 ;
            }


            }
            return -1 ;

        }

};

/*
general
1 1 1 0 0 0 -20 -20 -20 200	- A B C D E F G H I J
0 0 0 0 0 5	- cube reference point, length, width, height (0 indicates no clipping along this dimension)
0.0 0.0 1.0	- color
0.4 0.2 0.1 0.3	- ambient, diffuse, specular, recursive reflection coefficient
3		- shininess
*/

class General : public Object{

    public:
        double a,b,c,d,e,f,g,h,i,j ;

        General(double A,double B,double C,double D,double E,double F,
                double G,double H,double I,double J , Vector3D ref_point ,
                        double length_,double width_,double height_ ){
            //cout << "General Constructor "<< endl;
            a = A ; b = B ; c = C ; d = D ; e = E ; f = F ; g = G ; h = H ; i = I ; j= J ;
            reference_point = ref_point ;
            length = length_ ;
            width = width_ ;
            height = height_ ;
        }
        void draw(){
            //cout << "In General draw : " << endl;
            glColor3f(color[0],color[1],color[2]);

        }
        int clip(double t , Ray r)
        {
            Vector3D intersectionPoint ;
            intersectionPoint.x = r.start.x + t*r.dir.x ;
            intersectionPoint.y = r.start.y + t*r.dir.y ;
            intersectionPoint.z = r.start.z + t*r.dir.z ;

            //cout << intersectionPoint.z << endl;
            if(length != 0){
                if( reference_point.x > intersectionPoint.x ||
                ((reference_point.x + length) < intersectionPoint.x)){
                    return 0;
                }
            }
            if( width != 0){
                if( reference_point.y > intersectionPoint.y ||
                (reference_point.y + width) < intersectionPoint.y){
                    return 0;
                }
            }
            if( height != 0){
                if( reference_point.z > intersectionPoint.z ||
                (reference_point.z + height < intersectionPoint.z)){
                    //cout << intersectionPoint.z << endl;
                    return 0;
                }
            }
            return 1 ;
        }
        double intersect(Ray r, double *col, int level){

            //cout << "General interesction \n" << endl;
            double a_eqn = a*(r.dir.x)*(r.dir.x) + b*(r.dir.y)*(r.dir.y) + c*(r.dir.z)*(r.dir.z)+
                    d*(r.dir.x)*(r.dir.y) + e*(r.dir.x)*(r.dir.z) + f *(r.dir.y)*(r.dir.z) ;

            double b_eqn = 2*a*(r.start.x)*(r.dir.x) + 2*b*(r.start.y)*(r.dir.y) + 2*c*(r.start.z)*(r.dir.z)+
                              d*((r.start.y)*(r.dir.x) + (r.start.x)*(r.dir.y))
                            + e * ((r.start.x)*(r.dir.z) + (r.start.z)*(r.dir.x))
                            + f *((r.start.y)*(r.dir.z) + (r.start.z)*(r.dir.y))
                            + g * r.dir.x + h * r.dir.y + i * r.dir.z ;

            double c_eqn = a*(r.start.x)*(r.start.x) + b*(r.start.y)*(r.start.y) + c*(r.start.z)*(r.start.z)+
                        d*(r.start.x)*(r.start.y) + e*(r.start.x)*(r.start.z) + f*(r.start.y)*(r.start.z) +
                        g*(r.start.x) + h*(r.start.y) + i*(r.start.z) + j ;

            double d_square = b_eqn*b_eqn - 4*a_eqn*c_eqn ;
            if(d_square < 0){
                return -1 ;
            }
            double d = sqrt(d_square) ;
            double t1 = (-b_eqn+d)/(2*a_eqn);
            double t2 = (-b_eqn-d)/(2*a_eqn);

            if(clip(t1,r)==0)
            {
                t1 = -1 ;
            }

            if(clip(t2,r)==0)
            {
                t2 = -1 ;
            }
            double t ;
            if(t1< 0 && t2 <0){
                t = -1 ;
            }else if(t1<0 && t2>0){
                t = t2 ;
            }else if(t1>0 && t2<0 ){
                t = t1 ;
            }else if(t1>0 && t2>0 && t1>t2  ){
                t = t2 ;
            }else if(t1>0 && t2>0 && t1<t2 ){
                t = t1 ;
            }

            if(level == 0){
                return t ;
            }


            //phong

            Vector3D intersectionPoint ;
            intersectionPoint.x = r.start.x + t*r.dir.x ;
            intersectionPoint.y = r.start.y + t*r.dir.y ;
            intersectionPoint.z = r.start.z + t*r.dir.z ;

            double intersectionPointColor[3];
            intersectionPointColor[0] = color[0];
            intersectionPointColor[1] = color[1];
            intersectionPointColor[2] = color[2];

            col[0] = intersectionPointColor[0]*coEfficients[0] ;
            col[1] = intersectionPointColor[1]*coEfficients[0] ;
            col[2] = intersectionPointColor[2]*coEfficients[0] ;

            Vector3D normal ;
            normal.x = 2*a*intersectionPoint.x + d* intersectionPoint.y +
                        e*intersectionPoint.z + g ;
            normal.y = 2*b*intersectionPoint.y + d* intersectionPoint.x +
                        f*intersectionPoint.z + h ;
            normal.z = 2*c*intersectionPoint.z + e* intersectionPoint.x +
                        f*intersectionPoint.y + i ;
            normal.normalize();
            if(dot_product(normal,r.dir) < 0){
                normal = scaleVector(normal,-1);
            }

            for(int i = 0 ; i < lights.size() ; i++){
                Ray ray1(intersectionPoint,lights[i].light_pos);
                ray1.start = Vector3D(intersectionPoint.x+ sml*ray1.dir.x ,
                                  intersectionPoint.y+ sml*ray1.dir.y ,
                                  intersectionPoint.z+ sml*ray1.dir.z) ;//Vector3D ray1(subtract_vector(intersectionPoint,lights[i].light_pos));

                int obscured = 0 ;
                double tMin = -1 ;
                for(int k = 0 ; k < objects.size() ; k++){
                    t = objects[k]->intersect(ray1,color,0);
                    if(tMin < 0 && t > 0){
                        tMin = t ;
                    }else if(tMin > t && t > 0){
                        tMin = t ;
                    }
                }
                Vector3D distance_vec = subtract_vector(lights[i].light_pos, intersectionPoint);
                double distance = distance_vec.magnitude();

                if(tMin > 0 && distance > tMin){
                    obscured = 1 ;
                }

                //if(tMin > 0 && abs(t_this - tMin) < 0.00000001 )
                if(obscured == 0 )
                {
                    //cout << "Here" << endl;
                    double lambert = dot_product(normal,ray1.dir) > 0 ?
                                      dot_product(normal,ray1.dir) : 0;
                    //find reflected ray, rayr for rayl
                    //R = 2 (L.N)N – L
                    Ray ray_r ;
                    ray_r.start = intersectionPoint ;
                    //confused
                    ray_r.dir = subtract_vector(scaleVector(normal,2*lambert),ray1.dir);
                    //ray_r.dir = subtract_vector(ray1.dir,scaleVector(normal,2*lambert));
                    ray_r.dir.normalize() ;

                    // r k negative kore
                    Vector3D negative_r = scaleVector(r.dir,-1);
                    negative_r.normalize();
                    double phongValue =  dot_product(negative_r,ray_r.dir) > 0 ?
                                dot_product(negative_r,ray_r.dir) : 0;

                    //cout << "lambert : "<< lambert << "  phong : " << phongValue << endl;
                    col[0] += lights[i].color[0]*coEfficients[1]*lambert*intersectionPointColor[0] ;
                    col[0] += lights[i].color[0]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[0] ;

                    col[1] += lights[i].color[1]*coEfficients[1]*lambert*intersectionPointColor[1] ;
                    col[1] += lights[i].color[1]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[1] ;

                    col[2] += lights[i].color[2]*coEfficients[1]*lambert*intersectionPointColor[2] ;
                    col[2] += lights[i].color[2]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[2] ;

                    if(col[0] < 0){
                        col[0] = 0 ;
                    }
                    if(col[1] < 0){
                        col[1] = 0 ;
                    }
                    if(col[2] < 0){
                        col[2] = 0 ;
                    }
                    if(col[0] > 1){
                        col[0] = 1 ;
                    }
                    if(col[1] > 1){
                        col[1] = 1 ;
                    }
                    if(col[2] > 1){
                        col[2] = 1 ;
                    }

                }
            }

            if(level >= level_of_recursion ){
                //cout << "level 4 : "<< col[0] << " "<< col[1] << " " << col[2] << endl;
                return t ;
            }

            Ray reflected ;
            //R = 2 (L.N)N – L
            double scale = 2* dot_product(r.dir,normal);
            //reflected.dir = subtract_vector(scaleVector(normal,scale),r.dir);
            reflected.dir = subtract_vector(r.dir , scaleVector(normal,scale));

            reflected.dir.normalize();
            reflected.start = Vector3D(intersectionPoint.x + sml*reflected.dir.x ,
                              intersectionPoint.y+ sml*reflected.dir.y ,
                              intersectionPoint.z+ sml*reflected.dir.z);

            double t_ , tMin_ = -1 ;
            double *reflectedColor = new double[3];
            reflectedColor[0] = 0.0; reflectedColor[1] = 0.0 ; reflectedColor[2] = 0.0 ;
            Object *o ;
            for(int k = 0 ; k < objects.size() ; k++){
                t_ = objects[k]->intersect(reflected , reflectedColor,0);
                if(tMin_ < 0 && t_ > 0){
                    tMin_ = t_ ;
                    o = objects[k];
                }else if(tMin_ > t_ && t_ > 0){
                    tMin_ = t_ ;
                    o = objects[k];
                }
            }
            if(tMin_ > 0){
                o->intersect(reflected, reflectedColor ,level+1);
                //cout << "level : " << level << endl;

                //cout << reflectedColor[0] << " " << reflectedColor[1] << " " << reflectedColor[2] << endl;
                col[0] += reflectedColor[0] * coEfficients[3];
                col[1] += reflectedColor[1] * coEfficients[3];
                col[2] += reflectedColor[2] * coEfficients[3];
            }

            if(col[0] < 0){
                col[0] = 0 ;
            }
            if(col[1] < 0){
                col[1] = 0 ;
            }
            if(col[2] < 0){
                col[2] = 0 ;
            }
            if(col[0] > 1){
                col[0] = 1 ;
            }
            if(col[1] > 1){
                col[1] = 1 ;
            }
            if(col[2] > 1){
                col[2] = 1 ;
            }

        }
};


class Floor : public Object{
    public:
        double floorWidth = 1000;
        double tile_width = 20;

        Floor(double floor_Width, double tileWidth){
            //cout << "Floor Constructor "<< endl;
            Vector3D ref_point(-floorWidth/2,floorWidth/2,0) ;
            reference_point = ref_point ;
            length = tileWidth ;
            floorWidth = floor_Width;
            tile_width = tileWidth ;
        }

        void draw()
        {

            for(int i = -floorWidth/(tile_width*2) ; i < floorWidth/(tile_width*2) ; i++){
                for(int j = -floorWidth/(tile_width*2) ; j < floorWidth/(tile_width*2) ; j++){
                    if(((i)%2 == 0 && (j)%2 == 0) || ((i)%2 != 0 && (j)%2 != 0)){
                        glColor3f(0,0,0);
                        //if(i == 0 && j == 0 )
                            //cout << "black"<< endl;
                    }
                    else{
                        glColor3f(1,1,1);
                    }
                     glBegin(GL_POLYGON);
                          glVertex3f(i*tile_width, j*tile_width, 0.0);
                          glVertex3f((i+1)*tile_width, j*tile_width, 0.0);
                          glVertex3f((i+1)*tile_width, (j+1)*tile_width, 0.0);
                          glVertex3f(i*tile_width, (j+1)*tile_width, 0.0);
                    glEnd();
                }
            }
        }
        double intersect(Ray r, double *col, int level){
            Vector3D normal(0,0,1.0);
            double t = - dot_product(normal,r.start)/dot_product(normal,r.dir);
            if( t < 0){
                return -1 ;
            }
            Vector3D intersectionPoint ;
            intersectionPoint.x = r.start.x + t*r.dir.x ;
            intersectionPoint.y = r.start.y + t*r.dir.y ;
            intersectionPoint.z = r.start.z + t*r.dir.z ;
            if(intersectionPoint.x > -floorWidth/2 && intersectionPoint.x < floorWidth/2 &&
               intersectionPoint.y > -floorWidth/2 && intersectionPoint.y < floorWidth/2 )
            {
                int row = round((intersectionPoint.x - reference_point.x  )/20.0) ;
                int col_ = round((intersectionPoint.y - reference_point.y )/20.0) ;
                if((row%2 == 0 && col_%2 == 0) || (row%2 != 0 && col_%2 != 0) ){
                    color[0] = 0 ;
                    color[1] = 0 ;
                    color[2] = 0 ;
                }else{
                    color[0] = 1 ;
                    color[1] = 1 ;
                    color[2] = 1 ;
                }
                if(level == 0 )
                    return t ;


                Vector3D intersectionPoint ;
                intersectionPoint.x = r.start.x + t*r.dir.x ;
                intersectionPoint.y = r.start.y + t*r.dir.y ;
                intersectionPoint.z = r.start.z + t*r.dir.z ;

                double intersectionPointColor[3];
                intersectionPointColor[0] = color[0];
                intersectionPointColor[1] = color[1];
                intersectionPointColor[2] = color[2];

                col[0] = intersectionPointColor[0]*coEfficients[0] ;
                col[1] = intersectionPointColor[1]*coEfficients[0] ;
                col[2] = intersectionPointColor[2]*coEfficients[0] ;

                Vector3D normal(0.0,0.0,1.0);
                normal.normalize();

                for(int i = 0 ; i < lights.size() ; i++){
                    Ray ray1(intersectionPoint,lights[i].light_pos);
                    ray1.start = Vector3D(intersectionPoint.x+ sml*ray1.dir.x ,
                                          intersectionPoint.y+ sml*ray1.dir.y ,
                                          intersectionPoint.z+ sml*ray1.dir.z);


                    int obscured = 0 ;
                    double tMin = -1 ;
                    for(int k = 0 ; k < objects.size() ; k++){
                        t = objects[k]->intersect(ray1,color,0);
                        if(tMin < 0 && t > 0){
                            tMin = t ;
                        }else if(tMin > t && t > 0){
                            tMin = t ;
                        }
                    }
                    Vector3D distance_vec = subtract_vector(lights[i].light_pos, intersectionPoint);
                    double distance = distance_vec.magnitude();

                    if(tMin > 0 && distance > tMin){
                        obscured = 1 ;
                    }
                    if(obscured == 0 )
                    {
                        double lambert = dot_product(normal,ray1.dir) > 0 ?
                                          dot_product(normal,ray1.dir) : 0;
                        //find reflected ray, rayr for rayl
                        //R = 2 (L.N)N – L
                        Ray ray_r ;
                        ray_r.start = intersectionPoint ;
                        //confused
                        ray_r.dir = subtract_vector(scaleVector(normal,2*lambert),ray1.dir);
                        //ray_r.dir = subtract_vector(ray1.dir,scaleVector(normal,2*lambert));
                        ray_r.dir.normalize() ;

                        // r k negative kore
                        Vector3D negative_r = scaleVector(r.dir,-1);
                        negative_r.normalize();
                        double phongValue =  dot_product(negative_r,ray_r.dir) > 0 ?
                                    dot_product(negative_r,ray_r.dir) : 0;

                        //cout << "lambert : "<< lambert << "  phong : " << phongValue << endl;
                        col[0] += lights[i].color[0]*coEfficients[1]*lambert*intersectionPointColor[0] ;
                        col[0] += lights[i].color[0]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[0] ;

                        col[1] += lights[i].color[1]*coEfficients[1]*lambert*intersectionPointColor[1] ;
                        col[1] += lights[i].color[1]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[1] ;

                        col[2] += lights[i].color[2]*coEfficients[1]*lambert*intersectionPointColor[2] ;
                        col[2] += lights[i].color[2]*coEfficients[2]*pow(phongValue,shine)*intersectionPointColor[2] ;

                           if(col[0] < 0){
                            col[0] = 0 ;
                            }
                            if(col[1] < 0){
                                col[1] = 0 ;
                            }
                            if(col[2] < 0){
                                col[2] = 0 ;
                            }
                            if(col[0] > 1){
                                col[0] = 1 ;
                            }
                            if(col[1] > 1){
                                col[1] = 1 ;
                            }
                            if(col[2] > 1){
                                col[2] = 1 ;
                            }

                        }

                    }


                    if(level >= level_of_recursion)
                        return t ;

                    Ray reflected ;
                    double scale = 2* dot_product(r.dir,normal);
                    reflected.dir = subtract_vector(scaleVector(normal,scale),r.dir);
                    //reflected.dir = subtract_vector(r.dir , scaleVector(normal,scale));

                    reflected.dir.normalize();
                    reflected.start = Vector3D(intersectionPoint.x + sml*reflected.dir.x ,
                                      intersectionPoint.y+ sml*reflected.dir.y ,
                                      intersectionPoint.z+ sml*reflected.dir.z);

                    double t_ , tMin_ = -1 ;
                    double *reflectedColor = new double[3];
                    reflectedColor[0] = 0.0; reflectedColor[1] = 0.0 ; reflectedColor[2] = 0.0 ;
                    Object *o ;
                    for(int k = 0 ; k < objects.size() ; k++){
                        t_ = objects[k]->intersect(reflected,reflectedColor,0);
                        if(tMin_ < 0 && t_ > 0){
                            tMin_ = t_ ;
                            o = objects[k];
                        }else if(tMin_ > t_ && t_ > 0){
                            tMin_ = t_ ;
                            o = objects[k];
                        }
                    }
                    if(tMin_ > 0){
                        o->intersect(reflected, reflectedColor ,level+1);
                        //cout << reflectedColor[0] << " " << reflectedColor[1] << " " << reflectedColor[2] << endl;
                        col[0] += reflectedColor[0] * coEfficients[3];
                        col[1] += reflectedColor[1] * coEfficients[3];
                        col[2] += reflectedColor[2] * coEfficients[3];

                          if(col[0] < 0){
                            col[0] = 0 ;
                        }
                        if(col[1] < 0){
                            col[1] = 0 ;
                        }
                        if(col[2] < 0){
                            col[2] = 0 ;
                        }
                        if(col[0] > 1){
                            col[0] = 1 ;
                        }
                        if(col[1] > 1){
                            col[1] = 1 ;
                        }
                        if(col[2] > 1){
                            col[2] = 1 ;
                        }

                   }



            }

            //return -1 ;
        }
};

