using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace WOA
{
    public class Whale
    {
        public Vector coords;
        public double fitness = 0;
        public Whale() { coords = new Vector(0, 0); fitness = 0; }
        public Whale(double x, double y, double fitness)
        {
            coords = new Vector(x, y);
            this.fitness = fitness;
        }
    }
}
