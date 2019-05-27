import os
import numpy as np
import sys
from carboxAxis import calcMassCenter, calcNewCoor, newPosition

def readTypeInfo(file):
	sour = open(file, 'r')
	atomMasses = []
	atomPosition = []
	mass_dict = {}
	for line in sour:
		if line == 'Masses\n':
			while True:
				massLine = sour.readline()
				if massLine == '\n':
					if mass_dict == {}:
						continue
					else:
						break
				words = massLine.split()
				id = words[0]
				mass = words[1]
				mass_dict[id] = mass
		if line == 'Atoms\n':
			while True:
				typeLine = sour.readline()
				if typeLine == '\n':
					if atomMasses == []:
						continue
					else:
						break
				words = typeLine.split()
				type = words[2]
				position = (words[4],words[5],words[6])
				atomMasses.append(mass_dict[type])
				atomPosition.append(position)
	sour.close()
	return [tuple(map(float,p)) for p in atomPosition], list(map(float,atomMasses))


def calcPointCloud(file, featureDir):
	head, tail = os.path.split(file)
	id = tail[tail.find("linker") + len("linker") : tail.find(".lmpdat")]

	massToIdx = {1:0, 12:1, 14:2, 16:3}
	position_list, mass_list = readTypeInfo(file)

	massCenter = calcMassCenter(mass_list, position_list)
	transMatrix = calcNewCoor(position_list)
	for i in range(len(position_list)):
		position_list[i] = newPosition(position_list[i], transMatrix, massCenter)

	mass_list = list(map(lambda x: round(x), mass_list))
	position = np.array(position_list)
	mass = np.zeros(shape=(len(mass_list), 4))
	for i,x in enumerate(mass_list):
		mass[i][massToIdx[x]] = 1
	positionWithMass = np.concatenate((position, mass), axis=1).astype(np.float16)
	# padding = np.zeros(shape=(57-positionWithMass.shape[0], positionWithMass.shape[1]), dtype=np.float16)
	# positionWithMass = np.concatenate((positionWithMass, padding), axis=0)
	np.save(os.path.join(featureDir, id + ".npy"), positionWithMass)
	print(os.path.join(featureDir, id + ".npy") + " computed")


if __name__ == "__main__":
	file = sys.argv[1]
	featureDir = sys.argv[2]
	calcPointCloud(file, featureDir)


