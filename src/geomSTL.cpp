#include "geomSTL.h"
#include "definitions.h"
#define MYSCALE 0.25

Vec3::Vec3( void )
{
}

Vec3::Vec3( float X, float Y, float Z )
{
    this->X = X;
    this->Y = Y;
    this->Z = Z;
}

float Vec3::Length()
{
    return sqrt( X * X + Y * Y + Z * Z );
}

Vec3 Vec3::Vectors()
{
    Vec3 vec;
    float X = this->X;
    float Y = this->Y;
    float Z = this->Z;
    vec.X = X;
    vec.Y = Y;
    vec.Z = Z;
    return vec;
}

Vec3 Vec3::Normalize()
{
    Vec3 vector;
    float length = this->Length();

    if ( length != 0 )
    {
        vector.X = X / length;
        vector.Y = Y / length;
        vector.Z = Z / length;
    }

    return vector;
}

Vec3::~Vec3( void )
{
}

//***********************************************************************//
/*
    std::ostream& GeomSTL::operator<<( std::ostream &out, const point p )
    {
        out << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
        return out;
    }

    std::ostream& GeomSTL::operator<<( std::ostream &out, const triangle &t )
    {
        out << "---- TRIANGLE ----" << std::endl;
        out << t.normal << std::endl;
        out << t.v1 << std::endl;
        out << t.v2 << std::endl;
        out << t.v3 << std::endl;
        return out;
    }
*/
    float GeomSTL::parse_float( std::ifstream &s )
    {
        char f_buf[sizeof( float )];
        s.read( f_buf, 4 );
        float *fptr = (float *)f_buf;
        return *fptr;
    }

    point GeomSTL::parse_point( std::ifstream &s )
    {
        float x = parse_float( s );
        float y = parse_float( s );
        float z = parse_float( s );
        return point( x, y, z );
    }

    stl_data GeomSTL::parse_stl( const std::string &stl_path )
    {
        std::ifstream stl_file( stl_path.c_str(), std::ios::in | std::ios::binary );
        if ( !stl_file )
        {
            std::cout << "ERROR: COULD NOT READ FILE" << std::endl;
            assert( false );
        }

        char header_info[80] = "";
        char n_triangles[4];
        stl_file.read( header_info, 80 );
        stl_file.read( n_triangles, 4 );
        std::string h( header_info );
        stl_data info( h );
        unsigned int *r = (unsigned int *)n_triangles;
        unsigned int num_triangles = *r;
        for ( unsigned int i = 0; i < num_triangles; i++ )
        {
            auto normal = parse_point( stl_file );
            auto v1 = parse_point( stl_file );
            auto v2 = parse_point( stl_file );
            auto v3 = parse_point( stl_file );
            info.triangles.push_back( triangle( normal, v1, v2, v3 ) );
            char dummy[2];
            stl_file.read( dummy, 2 );
        }
        return info;
    }

//
//
//
void GeomSTL::construct(real* xyz1)
{
for(int i=0;i<3;i++)
{
xyz[i]=xyz1[i];
}

}



void GeomSTL::readSTLGeom( char *argv[], const real *xyz )
{
    real scale = 0.02;

    real disp = 1.0;

    auto info = parse_stl( argv[1] );
#if ( 1 )

    std::vector<triangle> triangles = info.triangles;
    // QMOut << "STL HEADER = " << info.name << std::endl;
    //  QMOut << "# triangles = " << triangles.size() << std::endl;

    real center[3] = {0.0, 0.0, 0.0};

    // This is for debugging
    // QMOut << info.triangles[88000]<<std::endl;
    // QMOut << triangles.size() << std::endl;

    // distribute the geometry between the blocks looking at their block extent (or ancestor length)

    for ( unsigned int i = 0; i < triangles.size(); i++ )
    {
        center[0] += ( info.triangles[i].v1.x + info.triangles[i].v2.x + info.triangles[i].v3.x );
        center[1] += ( info.triangles[i].v1.y + info.triangles[i].v2.y + info.triangles[i].v3.y );
        center[2] += ( info.triangles[i].v1.z + info.triangles[i].v2.z + info.triangles[i].v3.z );
    }

    center[0] = center[0] / triangles.size() / 3.0;
    center[1] = center[1] / triangles.size() / 3.0;
    center[2] = center[2] / triangles.size() / 3.0;

    // min and max can be found and then converted to center by division by 3

    real xmax = ( info.triangles[0].v1.x + info.triangles[0].v2.x + info.triangles[0].v3.x );
    real ymax = ( info.triangles[0].v1.y + info.triangles[0].v2.y + info.triangles[0].v3.y );
    real zmax = ( info.triangles[0].v1.z + info.triangles[0].v2.z + info.triangles[0].v3.z );

    for ( unsigned int i = 0; i < triangles.size(); i++ )
    {
        if ( info.triangles[i].v1.x + info.triangles[i].v2.x + info.triangles[i].v3.x > xmax )
        {
            xmax = info.triangles[i].v1.x + info.triangles[i].v2.x + info.triangles[i].v3.x;
        }

        if ( info.triangles[i].v1.y + info.triangles[i].v2.y + info.triangles[i].v3.y > ymax )
        {
            ymax = info.triangles[i].v1.y + info.triangles[i].v2.y + info.triangles[i].v3.y;
        }
        if ( info.triangles[i].v1.z + info.triangles[i].v2.z + info.triangles[i].v3.z > zmax )
        {
            zmax = info.triangles[i].v1.z + info.triangles[i].v2.z + info.triangles[i].v3.z;
        }
    }

    scale = xmax;
    if ( ymax > scale )
    {
        scale = ymax;
    }
    if ( zmax > scale )
    {
        scale = zmax;
    }
    scale = 3. / scale * MYSCALE;

    for ( unsigned int i = 0; i < triangles.size(); i++ )
    {
        info.triangles[i].v1.x = info.triangles[i].v1.x * scale - center[0] * scale;
        info.triangles[i].v1.y = info.triangles[i].v1.y * scale - center[1] * scale;
        info.triangles[i].v1.z = info.triangles[i].v1.z * scale - center[2] * scale;

        info.triangles[i].v2.x = info.triangles[i].v2.x * scale - center[0] * scale;
        info.triangles[i].v2.y = info.triangles[i].v2.y * scale - center[1] * scale;
        info.triangles[i].v2.z = info.triangles[i].v2.z * scale - center[2] * scale;

        info.triangles[i].v3.x = info.triangles[i].v3.x * scale - center[0] * scale;
        info.triangles[i].v3.y = info.triangles[i].v3.y * scale - center[1] * scale;
        info.triangles[i].v3.z = info.triangles[i].v3.z * scale - center[2] * scale;
    }

    geom_nn = 0;
    int bol[3] = {0, 0, 0};

    for ( unsigned int i = 0; i < triangles.size(); i++ )
    {
        bol[0] = ( info.triangles[i].v1.x + info.triangles[i].v2.x + info.triangles[i].v3.x ) <= xyz[1] * 3.
                 && ( info.triangles[i].v1.x + info.triangles[i].v2.x + info.triangles[i].v3.x ) >= xyz[0] * 3.;
        bol[1] = ( info.triangles[i].v1.y + info.triangles[i].v2.y + info.triangles[i].v3.y ) <= xyz[3] * 3.
                 && ( info.triangles[i].v1.y + info.triangles[i].v2.y + info.triangles[i].v3.y ) >= xyz[2] * 3.;
        bol[2] = ( info.triangles[i].v1.z + info.triangles[i].v2.z + info.triangles[i].v3.z ) <= xyz[5] * 3.
                 && ( info.triangles[i].v1.z + info.triangles[i].v2.z + info.triangles[i].v3.z ) >= xyz[4] * 3;

        if ( bol[0] && bol[1] && bol[2] )
        {
            ( geom_nn )++;
        }
    }
    //  printf( "nn=%d\n", geom_nn );

    triangle_center = (real *)malloc( 3 * ( geom_nn ) * sizeof( real ) );

    int count = 0;

    for ( unsigned int i = 0; i < triangles.size(); i++ )
    {
        bol[0] = ( info.triangles[i].v1.x + info.triangles[i].v2.x + info.triangles[i].v3.x ) <= xyz[1] * 3.
                 && ( info.triangles[i].v1.x + info.triangles[i].v2.x + info.triangles[i].v3.x ) >= xyz[0] * 3.;
        bol[1] = ( info.triangles[i].v1.y + info.triangles[i].v2.y + info.triangles[i].v3.y ) <= xyz[3] * 3.
                 && ( info.triangles[i].v1.y + info.triangles[i].v2.y + info.triangles[i].v3.y ) >= xyz[2] * 3.;
        bol[2] = ( info.triangles[i].v1.z + info.triangles[i].v2.z + info.triangles[i].v3.z ) <= xyz[5] * 3.
                 && ( info.triangles[i].v1.z + info.triangles[i].v2.z + info.triangles[i].v3.z ) >= xyz[4] * 3;

        if ( bol[0] && bol[1] && bol[2] )
        {
            ( triangle_center )[3 * count] = ( info.triangles[i].v1.x + info.triangles[i].v2.x + info.triangles[i].v3.x ) / 3.0;
            ( triangle_center )[3 * count + 1] = ( info.triangles[i].v1.y + info.triangles[i].v2.y + info.triangles[i].v3.y ) / 3.0;
            ( triangle_center )[3 * count + 2] = ( info.triangles[i].v1.z + info.triangles[i].v2.z + info.triangles[i].v3.z ) / 3.0;
            //	        cout<<(*triangle_center)[3*count ]<<"\t" <<(*triangle_center)[3*count+1]<<"\t"
            //<<(*triangle_center)[3*count+2]<<endl;
            count++;
        }
    }

    printf( "scale=%16.16lf %16.16lf %16.16lf %16.16lf \n", scale, center[0] * scale, center[1] * scale, center[2] * scale );

    if ( CHECK_MESH )
    {
        checkMesh( triangles );
    }

   geom_xyz=triangle_center;
#endif
}

void GeomSTL::checkMesh( std::vector<triangle> &triangles )
{
    // Area Normal Components
    float AreaX = 0.0f;
    float AreaY = 0.0f;
    float AreaZ = 0.0f;
    float AreaXSum = 0.0f;
    float AreaYSum = 0.0f;
    float AreaZSum = 0.0f;

    // Length of Sides
    float MaxEdge;
    float MinEdge;
    float GlobalMaxEdge = 0.0f;
    float GlobalMinEdge = 1.0f; // Initialize to Anything but zero

    // Angles
    float Angle1;
    float Angle2;
    float Angle3;
    float MaxAng;
    float MinAng;
    float pi = 3.1415926535897;
    float GlobalMaxAng = 0.0f;
    float GlobalMinAng = pi;

    // Aspect Ratio
    float s; // Semiperimeter
    float K; // Area K
    float r; // K/s
    float R; // (abc)/4K
    float AR;
    float avgAR = 0.0f;
    float JacX = 0.0f;
    float JacY = 0.0f;
    float JacZ = 0.0f;

    // Vectors to store Jacobian at each of the nodes
    float *Jac_Node1 = new float[3];
    float *Jac_Node2 = new float[3];
    float *Jac_Node3 = new float[3];

    // Dot Product at each vertex with Area Vector
    float Jac1_dot_Area = 0.0f;
    float Jac2_dot_Area = 0.0f;
    float Jac3_dot_Area = 0.0f;

    // Avg X, Y, Z components of Average Jacobian
    float Avg_JacX = 0.0f;
    float Avg_JacY = 0.0f;
    float Avg_JacZ = 0.0f;
    float AvgJac_dot_Area = 0.0f;

    // Keep track of the cells
    unsigned int Pos = 0;
    unsigned int PosSkew = 0;
    unsigned int Neg = 0;

    unsigned int NegSkew = 0;

    ofstream    QMOut;
    std::string filename = "STL_";
    filename.append( "Quality_Metrics" );
    QMOut.open( filename );
 
    QMOut << "\n \n *****************************************************************************************\n \n" << endl;
    QMOut << "                             STL Quality Metrics" << endl;
    QMOut << "\n \n *****************************************************************************************\n \n" << endl;

    // Loop to Cycle through all of the triangle elements
    for ( unsigned int i = 0; i < triangles.size(); i++ )
    {
        // Initialize Vectors for each of the 3 sides
        Vec3 v1( triangles[i].v2.x - triangles[i].v1.x, triangles[i].v2.y - triangles[i].v1.y, triangles[i].v2.z - triangles[i].v1.z );

        Vec3 v2( triangles[i].v3.x - triangles[i].v2.x, triangles[i].v3.y - triangles[i].v2.y, triangles[i].v3.z - triangles[i].v2.z );

        Vec3 v3( triangles[i].v1.x - triangles[i].v3.x, triangles[i].v1.y - triangles[i].v3.y, triangles[i].v1.z - triangles[i].v3.z );

        // Compute Area Vector
        AreaX = 0.5f * ( ( v1.Vectors().Y * v2.Vectors().Z ) - ( v1.Vectors().Z * v2.Vectors().Y ) );
        AreaY = 0.5f * ( ( v1.Vectors().Z * v2.Vectors().X ) - ( v1.Vectors().X * v2.Vectors().Z ) );
        AreaZ = 0.5f * ( ( v1.Vectors().X * v2.Vectors().Y ) - ( v1.Vectors().Y * v2.Vectors().X ) );

        /*      // ***********Debug****************
                if ( i == 0 ) {
                QMOut << "Node 1: (" << triangles[i].v1.x << ", " << triangles[i].v1.y << ", " << triangles[i].v1.z << ")" <<endl;

                QMOut << "Node 2: (" << triangles[i].v2.x << ", " << triangles[i].v2.y << ", " << triangles[i].v2.z << ")" <<endl;

                QMOut << "Node 3: (" << triangles[i].v3.x << ", " << triangles[i].v3.y << ", " << triangles[i].v3.z << ")" <<endl;

                QMOut << "VECTOR1: [" << v1.Vectors().X << ", " << v1.Vectors().Y << ", " << v1.Vectors().Z << "]" <<endl;

                QMOut << "VECTOR2: [" << v2.Vectors().X << ", " << v2.Vectors().Y << ", " << v2.Vectors().Z << "]" <<endl;

                QMOut << "VECTOR3: [" << v3.Vectors().X << ", " << v3.Vectors().Y << ", " << v3.Vectors().Z << "]" <<endl;


                QMOut << "Area Vector: [" << AreaX << ", " << AreaY<< ", " << AreaZ << "]"  << endl;
                        QMOut << " JAC 1: " << Jac_Node1[0] << endl;

                        QMOut << " JAC 2: " << Jac_Node1[1] << endl;

                        QMOut << " JAC 3: " << Jac_Node1[2] << endl;

        //              QMOut <<"Dot Product: " << Jac_dot_Area << endl;

                        QMOut << "Normal Vec 1: [" << v1.Normalize().X << ", " << v1.Normalize().Y << ", " << v1.Normalize().Z << "]" <<
        endl;

                        QMOut << "Normal Vec 2: [" << v2.Normalize().X << ", " << v2.Normalize().Y << ", " << v2.Normalize().Z << "]" <<
        endl;

                        QMOut << "Normal Vec 3: [" << v3.Normalize().X << ", " << v3.Normalize().Y << ", " << v3.Normalize().Z << "]" <<
        endl;
                }

                // ***************************************
        */
        // Aspect Ratio:
        s = 0.5 * ( v1.Length() + v2.Length() + v3.Length() );

        K = sqrt( s * ( s - v1.Length() ) * ( s - v2.Length() ) * ( s - v3.Length() ) );

        r = K / s;

        R = ( v1.Length() * v2.Length() * v3.Length() ) / ( 4.0f * K );

        AR = R / ( 2.0f * r );

        // Length Skewness
        MaxEdge = std::max( v1.Length(), std::max( v2.Length(), v3.Length() ) );
        MinEdge = std::min( v1.Length(), std::min( v2.Length(), v3.Length() ) );

        GlobalMaxEdge = std::max( MaxEdge, GlobalMaxEdge );
        GlobalMinEdge = std::min( MinEdge, GlobalMinEdge );

        // Angle Skewness
        // Law of Cosines
        Angle3 = acos( ( ( v1.Length() * v1.Length() ) + ( v2.Length() * v2.Length() ) - ( v3.Length() * v3.Length() ) )
                       / ( 2.0f * v1.Length() * v2.Length() ) );

        // Law of Sines
        Angle2 = asin( ( v2.Length() * sin( Angle3 ) ) / ( v3.Length() ) );

        // Finding the Last Angle
        Angle1 = pi - Angle2 - Angle3;

        // Compare and save largest and smallest angles
        MaxAng = std::max( Angle1, std::max( Angle2, Angle3 ) );
        MinAng = std::min( Angle1, std::max( Angle2, Angle3 ) );

        GlobalMaxAng = std::max( MaxAng, GlobalMaxAng );
        GlobalMinAng = std::min( MinAng, GlobalMinAng );

        if ( GlobalMinAng == 0 )
        {
            QMOut << "Angle = 0" << endl;
            QMOut << "Element #: " << i << endl;
        }
        // Node 1
        // X-Component
        Jac_Node1[0] = ( v1.Normalize().Y * ( -v3.Normalize().Z ) ) - ( v1.Normalize().Z * ( -v2.Normalize().Y ) );

        // Y-Component
        Jac_Node1[1] = ( v1.Normalize().Z * ( -v3.Normalize().X ) ) - ( v1.Normalize().X * ( -v3.Normalize().Z ) );

        // Z-Component
        Jac_Node1[2] = ( v1.Normalize().X * ( -v3.Normalize().Y ) ) - ( v1.Normalize().Y * ( -v3.Normalize().X ) );

        // Dot Product with Normal
        Jac1_dot_Area = ( Jac_Node1[0] * AreaX ) + ( Jac_Node1[1] * AreaY ) + ( Jac_Node1[2] * AreaZ );

        // Node 2
        // X-Component
        Jac_Node2[0] = ( v2.Normalize().Y * ( -v1.Normalize().Z ) ) - ( v2.Normalize().Z * ( -v1.Normalize().Y ) );

        // Y-Component
        Jac_Node2[1] = ( v2.Normalize().Z * ( -v1.Normalize().X ) ) - ( v2.Normalize().X * ( -v1.Normalize().Z ) );

        // Z-Component
        Jac_Node2[2] = ( v2.Normalize().X * ( -v1.Normalize().Y ) ) - ( v2.Normalize().Y * ( -v1.Normalize().X ) );

        // Dot Product with Normal
        Jac2_dot_Area = ( Jac_Node2[0] * AreaX ) + ( Jac_Node2[1] * AreaY ) + ( Jac_Node2[2] * AreaZ );

        // Node 3
        // X-Component
        Jac_Node3[0] = ( v3.Normalize().Y * ( -v2.Normalize().Z ) ) - ( v3.Normalize().Z * ( -v2.Normalize().Y ) );

        // Y-Component
        Jac_Node3[1] = ( v3.Normalize().Z * ( -v2.Normalize().X ) ) - ( v3.Normalize().X * ( -v2.Normalize().Z ) );

        // Z-Component
        Jac_Node3[2] = ( v3.Normalize().X * ( -v2.Normalize().Y ) ) - ( v3.Normalize().Y * ( -v2.Normalize().X ) );

        // Dot Product with Normal
        Jac3_dot_Area = ( Jac_Node3[0] * AreaX ) + ( Jac_Node3[1] * AreaY ) + ( Jac_Node3[2] * AreaZ );

        Avg_JacX = ( Jac_Node1[0] + Jac_Node2[0] + Jac_Node3[0] ) / 3.0f;
        Avg_JacY = ( Jac_Node1[1] + Jac_Node2[1] + Jac_Node3[1] ) / 3.0f;
        Avg_JacZ = ( Jac_Node1[2] + Jac_Node2[2] + Jac_Node3[2] ) / 3.0f;

        AvgJac_dot_Area = ( Avg_JacX * AreaX ) + ( Avg_JacY * AreaY ) + ( Avg_JacZ * AreaZ );

        // Check for really deformed elements (porsche has one)
        if ( ( isinf( AR ) ) || ( isnan( AR ) ) )
        {
            QMOut << ">>>>>>>>>>>>>>>>>>>> Warning!!! There is an element with negative Aspect Ratio <<<<<<<<<<<<<<<<<<<<" << endl;
            QMOut << "Element #: " << i << "\t"
                      << "avgAR: " << avgAR << endl;
            QMOut << "s: " << s << endl;
            QMOut << "K: " << K << endl;
            QMOut << "R: " << R << endl;
            QMOut << "AR: " << AR << endl;
            QMOut << "Vector 1: " << v1.Length() << endl;
            QMOut << "Vector 2: " << v2.Length() << endl;
            QMOut << "Vector 3: " << v3.Length() << endl;
            QMOut << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
            continue;
        }

        // Sum
        avgAR = avgAR + AR;

        AreaXSum += AreaX;
        AreaYSum += AreaY;
        AreaZSum += AreaZ;
        if ( ( Jac1_dot_Area > 0 ) && ( Jac2_dot_Area > 0 ) && ( Jac3_dot_Area > 0 ) )
        {
            // QMOut << "The cell is tagged as Positive" << endl;
            Pos += 1.0f;
        }
        else if ( ( Jac1_dot_Area < 0 ) && ( Jac2_dot_Area < 0 ) && ( Jac3_dot_Area < 0 ) )
        {
            // QMOut << "The cell is tagged as Negative" << endl;
            Neg += 1.0f;
        }
        else if ( ( ( Jac1_dot_Area < 0 ) && ( AvgJac_dot_Area > 0 ) ) || ( ( Jac2_dot_Area < 0 ) && ( AvgJac_dot_Area > 0 ) )
                  || ( ( Jac3_dot_Area < 0 ) && ( AvgJac_dot_Area > 0 ) ) )
        {
            // QMOut << "The cell is tagged as Positive Skewed" << endl;
            PosSkew += 1.0f;
        }
        else if ( ( ( Jac1_dot_Area > 0 ) && ( AvgJac_dot_Area < 0 ) ) || ( ( Jac2_dot_Area > 0 ) && ( AvgJac_dot_Area < 0 ) )
                  || ( ( Jac3_dot_Area > 0 ) && ( AvgJac_dot_Area < 0 ) ) )
        {
            // QMOut << "The cell is tagged as Negative Skewed\n" << endl;
            NegSkew += 1.0f;
        }
    }
    // Print Stats
    QMOut << "\nAverage Aspect Ratio: " << avgAR / triangles.size() << endl;
    QMOut << "\n\nLength Skewness: " << GlobalMaxEdge / GlobalMinEdge << endl;

    QMOut << "\n\tMax Edge: " << GlobalMaxEdge << "\t"
              << "Min Edge: " << GlobalMinEdge << endl;

    QMOut << "\n\nAngle Skewness: " << GlobalMaxAng / GlobalMinAng << endl;

    QMOut << "\n\tMax Angle (degrees): " << GlobalMaxAng *( 180.0f / pi ) << "\t"
              << "Min Angle (degrees): " << GlobalMinAng * ( 180.0f / pi ) << endl;
    QMOut << "\n\nArea Vector Sum: [" << AreaXSum << ", " << AreaYSum << ", " << AreaZSum << "]" << endl;

    QMOut << "\n\nNumber of Positive Cells (All Jacobians are Positive): " << Pos << endl;
    QMOut << "\nNumber of Positive Skewed Cells (Average of Jacobians is Positive): " << PosSkew << endl;
    QMOut << "\nNumber of Negative Cells (All Jacobians are Negative): " << Neg << endl;
    QMOut << "\nNumber of Negative Skewed Cells (Average of Jacobians is Negative): " << NegSkew << endl;
    QMOut << "\n******************************************************************************************" << endl;

    QMOut.close();
}
