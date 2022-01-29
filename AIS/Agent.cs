using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace DA
{
    public class Agent
    {
        public Vector coords;
        public double fitness = 0;
        public Agent() { coords = new Vector(0, 0); fitness = 0; }
        public Agent(double x, double y, double fitness)
        {
            coords = new Vector(x, y);
            this.fitness = fitness;
        }
    }
}
