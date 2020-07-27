using System;
using System.Linq;
using System.Xml.Schema;
using GraphSynth.Representation;
using System.Collections.Generic;
using OpenBabelFunctions;
using OpenBabel;
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

        public static double[] FindLengthAndRadius(designGraph host)
        {
            var LengthAndRadius = new double[2];
            var BeginningAtom = new double[3];
            var EndingAtom = new double[3];
            int m = 0;
            foreach (var n in host.nodes)
            {
                if (!n.localLabels.Contains("Boundary")) continue; 
                m += 1;
                if (m == 1)
                {
                    BeginningAtom[0] = n.X;
                    BeginningAtom[1] = n.Y;
                    BeginningAtom[2] = n.Z;
                }
                if (m == 2)
                {
                    EndingAtom[0] = n.X;
                    EndingAtom[1] = n.Y;
                    EndingAtom[2] = n.Z;
                }
            }
            //Length of the cylinder
            //Math.Pow((BeginningAtom[0] - EndingAtom[0]),2)
            LengthAndRadius[0] = Math.Sqrt(Math.Pow((BeginningAtom[0] - EndingAtom[0]), 2) + Math.Pow((BeginningAtom[1] - EndingAtom[1]), 2) + Math.Pow((BeginningAtom[2] - EndingAtom[2]), 2));

            //var d = BeginningAtom.subtract(EndingAtom);
            //LengthAndRadius[0]=d.norm2();

            var AxisVector = new double[3];
            AxisVector[0] = -(BeginningAtom[0] - EndingAtom[0]);
            AxisVector[1] = -(BeginningAtom[1] - EndingAtom[1]);
            AxisVector[2] = -(BeginningAtom[2] - EndingAtom[2]);


            var max = new double();
            foreach (var n in host.nodes)
            {
                var FarestPointToBeginningVector = new double[3];
                var area = new double();


                FarestPointToBeginningVector[0] = n.X - BeginningAtom[0];
                FarestPointToBeginningVector[1] = n.Y - BeginningAtom[1];
                FarestPointToBeginningVector[2] = n.Z - BeginningAtom[2];

                //(a1,a2,a3)x(b1,b2,b3)=(a2b3-a3b2,a3b1-a1b3,a1b2-a2b1)
                // FarestPointToBeginningVector[1]* AxisVector[2]- FarestPointToBeginningVector[2]* AxisVector[1]
                //FarestPointToBeginningVector[2] * AxisVector[0] - FarestPointToBeginningVector[0] * AxisVector[2]
                //FarestPointToBeginningVector[0] * AxisVector[1] - FarestPointToBeginningVector[1] * AxisVector[0]
                area = Math.Sqrt(Math.Pow(FarestPointToBeginningVector[1] * AxisVector[2] - FarestPointToBeginningVector[2] * AxisVector[1], 2) + Math.Pow(FarestPointToBeginningVector[2] * AxisVector[0] - FarestPointToBeginningVector[0] * AxisVector[2], 2) + Math.Pow(FarestPointToBeginningVector[0] * AxisVector[1] - FarestPointToBeginningVector[1] * AxisVector[0], 2));
                max = area / LengthAndRadius[0];
                if (LengthAndRadius[1] < max)
                    LengthAndRadius[1] = max;

            }

            return LengthAndRadius;
        }

        public static double[] CalcMoment(candidate cand)
        {
            var result = new double[10];

            // call openBabel functions to get x, y, z's into
            var mol = OBFunctions.designgraphtomol(cand.graph);
            //var newMol = OBFunctions.InterStepMinimize(mol);
            //OBFunctions.updatepositions(cand.graph, newMol);

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