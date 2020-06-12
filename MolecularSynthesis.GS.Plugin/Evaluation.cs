using System;
using System.Linq;
using System.Xml.Schema;
using GraphSynth.Representation;
using System.Collections.Generic;
//using GraphMolWrap;


namespace MolecularSynthesis.Plugin
{
    public static class Evaluation
    {
        static readonly Dictionary<string, double> massTabel = new Dictionary<string, double>()
        {
            {"H", 1},
            {"C", 12},
            {"N", 14},
            {"O", 16}
        };
        public static int CountAtoms(candidate cand)
        {
            return cand.graph.nodes.Count;
        }

        public static double[] CalcMoment(candidate cand)
        {
            var result = new double[10];
            // call openBabel functions to get x, y, z's into

            foreach (var n in cand.graph.nodes)
            {
                var mass = GetAtomMass(n);// OBFunctions.
                result[0] += mass;

                result[1] += mass * n.X;
                result[2] += mass * n.Y;
                result[3] += mass * n.Z;

                result[4] += mass * n.X * n.X;
                result[5] += mass * n.Y * n.Y;
                result[6] += mass * n.Z * n.Z;
                result[7] += mass * n.X * n.Y;
                result[8] += mass * n.Y * n.Z;
                result[9] += mass * n.Z * n.X;
            }

            return result;
        }

        private static double GetAtomMass(node n)
        {
            foreach (var label in n.localLabels)
            {
                if (massTabel.ContainsKey(label))
                    return massTabel[label];
            }

            throw new Exception("Error!: Node does not have a element label!");
        }

        internal static double distance(candidate child, double[] desiredMoment)
        {
            var childMoment = CalcMoment(child);
            //var difference = StarMath.norm1(StarMath.subtract(desiredMoment, childMoment));
            // this can be simplified using extensions...
            //var difference = desiredMoment.subtract(childMoment).norm1()/desiredMoment.add(childMoment).norm1();
            double[] difference = new double[10];
            for (int i = 0; i < 10; i++)
            {
                difference[i] = Math.Abs(desiredMoment[i] - childMoment[i]) / Math.Abs(desiredMoment[i] + childMoment[i]);
            }

            return difference.norm1(); 
        }
        public static double norm1(this IEnumerable<double> x)
        {
            return x.Sum(Math.Abs);
        }

    }
}