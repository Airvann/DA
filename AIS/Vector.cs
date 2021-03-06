using System;

namespace DA
{
    public class Vector
    {
        public double[] vector = new double[2];
        public static int dim = 2;

        public Vector(double x, double y)
        {
            vector[0] = x;
            vector[1] = y;
        }
        public Vector() { }

        public static Vector Norm(Vector vector)
        {
            Vector tmp = new Vector();
            for (int i = 0; i < dim; i++)
                tmp[i] = Math.Abs(vector[i]);
            return tmp;
        }

        public static Vector operator -(Vector vec1, Vector vec2)
        {
            Vector tmp = new Vector();
            tmp[0] = vec1[0] - vec2[0];
            tmp[1] = vec1[1] - vec2[1];
            return tmp;
        }

        public static Vector operator -(Vector vec)
        {
            Vector tmp = new Vector();
            tmp[0] = -vec[0];
            tmp[1] = -vec[1];
            return tmp;
        }

        public static Vector operator +(Vector vec1, Vector vec2)
        {
            Vector tmp = new Vector();
            tmp[0] = vec1[0] + vec2[0];
            tmp[1] = vec1[1] + vec2[1];
            return tmp;
        }
        public static Vector operator *(Vector vec1, double val)
        {
            Vector tmp = new Vector();
            tmp[0] = vec1[0] * val;
            tmp[1] = vec1[1] * val;
            return tmp;
        }

        public static Vector operator *(Vector vec1, Vector vec2)
        {
            Vector tmp = new Vector();
            tmp[0] = vec1[0] * vec2[0];
            tmp[1] = vec1[1] * vec2[0];
            return tmp;
        }

        public static Vector operator *(double val, Vector vec1)
        {
            return vec1 * val;
        }
        public static Vector operator +(Vector vec1, double val)
        {
            Vector tmp = new Vector();
            tmp[0] = vec1[0] + val;
            tmp[1] = vec1[1] + val;
            return tmp;
        }

        public static Vector operator +(double val, Vector vec1)
        {
            return vec1 + val;
        }
        public static Vector operator /(Vector vec1, double val)
        {
            Vector tmp = new Vector();
            tmp[0] = vec1[0] / val;
            tmp[1] = vec1[1] / val;
            return tmp;
        }

        public double this[int index]
        {
            get{ return vector[index]; }
            set{ vector[index] = value; }
        }
    }
}
