using OpenBabel;
using System;
using OpenBabelFunctions;
using GraphSynth.Search;
using GraphSynth;
using MolecularSynthesis.Plugin;
using System.Xml.XPath;
using System.IO;


namespace TestOpenBabel
{
    public class TestBenzRing : SearchProcess
    {
        public override string text => "TestBenz";
        public string filename = @"..\..\..\..\test";
        public string extension = "mol";
        //deault constructor
        public TestBenzRing(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }
        protected override void Run()
        {
            var resultMol = OBFunctions.designgraphtomol(seedGraph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(seedGraph, resultMol);
            SearchIO.addAndShowGraphWindow(seedGraph);

            var result = new double[2];
            result=Evaluation.FindLengthAndRadius(seedGraph);
            //[0] is Length, [1] is Radius
            SearchIO.output(result[0]+ " " +result[1]);

            var conv = new OBConversion();
            conv.SetInAndOutFormats("pdb", extension);
            conv.WriteFile(resultMol, filename + "." + extension);
            //File.AppendText()
        }
    }




}

