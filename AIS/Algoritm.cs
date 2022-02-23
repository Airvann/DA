using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics;

namespace DA
{
    //Класс алгоритма
    public class Algoritm
    {
        //Генератор псевдослучайных чисел
        private Random rand = new Random();

        #region Параметры алгоритма
        //Размер популяции стрекоз
        public int population;

        //Номер выбранной функции
        public int F;

        //область определения
        public double[,] D;

        //Максимальное число итераций
        public int MaxCount { get; set; }

        //разделение стрекоз в стае
        public double s;

        //выравнивание стрекоз в стае
        public double a;

        //сплоченность стрекоз в стае
        public double c;

        //стремление к лучшему решению
        public double f;

        //уклонение от худшего решения
        public double e;

        //память о предыстории
        public double w;

        //Начальное значение радиуса окрестности
        public Vector R = new Vector();
        //Текущая итерация
        public int currentIteration = 0;

        public Vector lb;
        public Vector ub;

        public Vector MaxDelta;
        #endregion

        //наиболее приспособленные особи
        public Agent best = new Agent();
        public Agent worst = new Agent();
        //Массив средней приспособленности
        public List<double> averageFitness = new List<double>();
        //Массив лучшей приспособленности
        public List<double> bestFitness = new List<double>();
        //Популяция стрекоз
        public List<Agent> individuals = new List<Agent>();
        //Pool
        public List<Agent> pool = new List<Agent>();
        //Шаги
        public List<Vector> steps;
        //Конструктор по умолчанию
        public Algoritm(){}
        //Старт алгоритма
        public Agent FastStartAlg(int population, int MaxCount, double[,] D, int F)
        {
            this.MaxCount = MaxCount;
            this.population = population;
            this.D = D;
            this.F = F;
            lb = new Vector(D[0, 0], D[1, 0]);
            ub = new Vector(D[0, 1], D[1, 1]);
            MaxDelta = (ub - lb) / 10;

            FormingPopulation();
            SetZeros();
            currentIteration = 0;

            for (; currentIteration < MaxCount; currentIteration++)
            {
                UpdateParams(currentIteration);
                PopulationOrder();
                NewPackGeneration();
                currentIteration++;
            }
            return PoolBest();
        }

        public void UpdateParams(int k)
        {
            w = 0.9 - k * ((0.9 - 0.4) / MaxCount);
            double my_c = 0.1 - k * ((0.9 - 0.4) / (0.5 * MaxCount));
            if (my_c < 0)
                my_c = 0;
            s = 2 * rand.NextDouble() * my_c;
            a = 2 * rand.NextDouble() * my_c;
            c = 2 * rand.NextDouble() * my_c;
            f = 2 * rand.NextDouble();
            e = my_c;
        }

        //Начальное формирование популяции
        public void FormingPopulation()
        {
            for (int i = 0; i < population; i++)
            {
                double x = rand.NextDouble();
                double y = rand.NextDouble();

                x = ((Math.Abs(D[0, 0]) + Math.Abs(D[0, 1])) * x - Math.Abs(D[0, 0]));
                y = ((Math.Abs(D[1, 0]) + Math.Abs(D[1, 1])) * y - Math.Abs(D[1, 0]));

                Agent Agent = new Agent(x, y, function(x, y, F));
                individuals.Add(Agent);
            }
        }

        public void PopulationOrder()
        {
            individuals = individuals.OrderBy(s => s.fitness).ToList();
            best = new Agent(individuals[0].coords[0], individuals[0].coords[1], individuals[0].fitness);
            worst = new Agent(individuals[population - 1].coords[0], individuals[population - 1].coords[1], individuals[population - 1].fitness);
            pool.Add(best);
        }

        //Формирование новой стаи
        public void NewPackGeneration()
        {
            R = ((ub - lb) / 4) + ((ub - lb) * (currentIteration / (0.5*MaxCount)));

            for (int i = 0; i < population; i++)
            {              
                List<Agent> neighbourhood = new List<Agent>();
                
                for (int j = 0; j < population; j++)
                {
                    double dist1 = Math.Abs(individuals[i].coords[0] - individuals[j].coords[0]);
                    double dist2 = Math.Abs(individuals[i].coords[1] - individuals[j].coords[1]);

                    if (dist1 < R[0] && dist2 < R[1] && dist1 != 0 && dist2 != 0)
                        neighbourhood.Add(individuals[j]);
                }
                if (neighbourhood.Count > 1)
                {
                    //Разделение
                    Vector S = new Vector(0, 0);
                    for (int m = 0; m < neighbourhood.Count; m++)
                        S += -(neighbourhood[m].coords - individuals[i].coords);

                    //Выравнивание
                    Vector A = new Vector(0, 0);
                    for (int m = 0; m < neighbourhood.Count; m++)
                        A += steps[m];
                    A /= neighbourhood.Count;

                    //Сплоченность
                    Vector C = new Vector(0, 0);
                    for (int m = 0; m < neighbourhood.Count; m++)
                        C += neighbourhood[m].coords - individuals[i].coords;
                    C /= neighbourhood.Count;

                    Vector F = best.coords - individuals[i].coords;
                    Vector E = worst.coords + individuals[i].coords;

                    steps[i] = s * S + a * A + c * C + f * F + e * E + w * steps[i];

                    if (steps[i][0] > MaxDelta[0])
                        steps[i][0] = MaxDelta[0];

                    if (steps[i][0] < -MaxDelta[0])
                        steps[i][0] = -MaxDelta[0];

                    if (steps[i][1] > MaxDelta[1])
                        steps[i][1] = MaxDelta[1];

                    if (steps[i][1] < -MaxDelta[1])
                        steps[i][1] = -MaxDelta[1];

                    individuals[i].coords += steps[i];

                    double x = individuals[i].coords[0];
                    double y = individuals[i].coords[1];

                    if (x < D[0, 0])
                        individuals[i].coords[0] = D[0, 0];
                    if (x > D[0, 1])
                        individuals[i].coords[0] = D[0, 1];
                    if (y < D[1, 0])
                        individuals[i].coords[1] = D[1, 0];
                    if (y > D[1, 1])
                        individuals[i].coords[1] = D[1, 1];
                }
                else 
                {
                    double beta = 1.5f;
                    double sigma = Math.Pow(SpecialFunctions.Gamma(1+beta)*Math.Sin(Math.PI*beta/2f)/(SpecialFunctions.Gamma((1+beta)/2f)*beta*Math.Pow(2,(beta-1)/2f)), 1/beta);
                    int dim = 2;
                    Vector Levy = new Vector();
                    Levy.vector[0] = 0.01 * rand.NextDouble() * dim * sigma / Math.Pow(Math.Abs(rand.NextDouble() * dim), 1 / beta);
                    Levy.vector[1] = 0.01 * rand.NextDouble() * dim * sigma / Math.Pow(Math.Abs(rand.NextDouble() * dim), 1 / beta);
                    individuals[i].coords = individuals[i].coords + Levy * individuals[i].coords;
                    SetZeros();
                }
                individuals[i].fitness = function(individuals[i].coords[0], individuals[i].coords[1], F);
            }
        }

        public void SetZeros()
        {
            steps = new List<Vector>();
            for (int i = 0; i < population; i++)
            {
                steps.Add(new Vector());
                steps[i][0] = 0;
                steps[i][1] = 0;
            }
        }
        private Agent PoolBest() 
        {
            pool = pool.OrderBy(s => s.fitness).ToList();
            return pool[0];
        }   
        //Все тестовые функции
        private float function(double x1, double x2, int F)
        {
            float funct = 0;
            if (F == 0)
                funct = -(float)(x1 * Math.Sin(Math.Sqrt(Math.Abs(x1))) + x2 * Math.Sin(Math.Sqrt(Math.Abs(x2))));
            else if (F == 1)
                funct = -(float)(x1 * Math.Sin(4 * Math.PI * x1) - x2 * Math.Sin(4 * Math.PI * x2 + Math.PI) + 1);
            else if (F == 2)
            {
                double[] c6 = Cpow(x1, x2, 6);
                funct = -(float)(1 / (1 + Math.Sqrt((c6[0] - 1) * (c6[0] - 1) + c6[1] * c6[1])));
            }
            else if (F == 3)
                funct = -(float)(0.5 - (Math.Pow(Math.Sin(Math.Sqrt(x1 * x1 + x2 * x2)), 2) - 0.5) / (1 + 0.001 * (x1 * x1 + x2 * x2)));
            else if (F == 4)
                funct = -(float)((-x1 * x1 + 10 * Math.Cos(2 * Math.PI * x1)) + (-x2 * x2 + 10 * Math.Cos(2 * Math.PI * x2)));
            else if (F == 5)
                funct = -(float)(-Math.E + 20 * Math.Exp(-0.2 * Math.Sqrt((x1 * x1 + x2 * x2) / 2)) + Math.Exp((Math.Cos(2 * Math.PI * x1) + Math.Cos(2 * Math.PI * x2)) / 2));
            else if (F == 6)
                funct = -(float)(Math.Pow(Math.Cos(2 * x1 * x1) - 1.1, 2) + Math.Pow(Math.Sin(0.5 * x1) - 1.2, 2) - Math.Pow(Math.Cos(2 * x2 * x2) - 1.1, 2) + Math.Pow(Math.Sin(0.5 * x2) - 1.2, 2));
            else if (F == 7)
                funct = -(float)(-Math.Sqrt(Math.Abs(Math.Sin(Math.Sin(Math.Sqrt(Math.Abs(Math.Sin(x1 - 1))) + Math.Sqrt(Math.Abs(Math.Sin(x2 + 2))))))) + 1);
            else if (F == 8)
                funct = -(float)(-(1 - x1) * (1 - x1) - 100 * (x2 - x1 * x1) * (x2 - x1 * x1));
            else if (F == 9)
                funct = -(float)(-x1 * x1 - x2 * x2);

            return funct;
        }
        //Вспомогательная функция
        private double[] Cpow(double x, double y, int p)
        {
            double[] Cp = new double[2];
            Cp[0] = x; Cp[1] = y;
            double x0 = 0;
            double y0 = 0;
            for (int i = 1; i < p; i++)
            {
                x0 = Cp[0] * x - Cp[1] * y;
                y0 = Cp[1] * x + Cp[0] * y;
                Cp[0] = x0; Cp[1] = y0;
            }
            return Cp;
        }
        //метод вычисления средней приспособленности
        public double AverageFitness()
        {
            double sum = 0;
            for (int i = 0; i < population; i++)
                sum += individuals[i].fitness;
            double fitness = (sum / population);
            averageFitness.Add(fitness);
            return fitness;
        }
    }
}
