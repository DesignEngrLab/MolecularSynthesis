import torch
import sys
from torch.utils.data import Dataset

if torch.cuda.is_available():
	sys.stdout.writelines("Using GPU")
	Tensor = torch.cuda.FloatTensor
else:
	sys.stdout.writelines("Using CPU")
	Tensor = torch.FloatTensor


class Tnet(torch.nn.Module):
	def __init__(self, K):
		super(Tnet, self).__init__()
		self.K = K
		self.fc1 = torch.nn.Linear(in_features=K, out_features=64)
		self.fc2 = torch.nn.Linear(in_features=64, out_features=128)
		self.fc3 = torch.nn.Linear(in_features=128, out_features=256)

		self.fc4 = torch.nn.Linear(in_features=256, out_features=256)
		self.fc5 = torch.nn.Linear(in_features=256, out_features=64)
		self.fc6 = torch.nn.Linear(in_features=64, out_features=K*K)

	def forward(self, x):
		self.max_pool_1d = torch.nn.MaxPool1d(kernel_size=x.shape[0])
		x = torch.relu(self.fc1(x))
		x = torch.relu(self.fc2(x))
		x = torch.relu(self.fc3(x))
		x, _ = torch.max(x, dim=0)
		x = torch.relu(self.fc4(x))
		x = torch.relu(self.fc5(x))
		x = self.fc6(x)
		x = x.reshape(self.K, self.K)
		return x



class PointNet(torch.nn.Module):
	def __init__(self, task):
		super(PointNet, self).__init__()
		self.transform1 = Tnet(3)
		self.fc1 = torch.nn.Linear(in_features=3, out_features=32)
		self.fc2 = torch.nn.Linear(in_features=32, out_features=64)
		self.transform2 = Tnet(68)
		self.fc3 = torch.nn.Linear(in_features=68, out_features=128)
		self.fc4 = torch.nn.Linear(in_features=128, out_features=256)

		self.fc5 = torch.nn.Linear(in_features=256, out_features=256)
		self.fc6 = torch.nn.Linear(in_features=256, out_features=64)
		if task == "Classification":
			self.fc7 = torch.nn.Linear(in_features=64, out_features=2)
			self.criterion = torch.nn.BCEWithLogitsLoss()
		elif task == "Regression":
			self.fc7 = torch.nn.Linear(in_features=64, out_features=1)
			self.criterion = torch.nn.MSELoss()
			
				
		self.optimizer = torch.optim.Adam(params=self.parameters())
		if torch.cuda.is_available():
			self.cuda()
		
		
		
	def forward(self, batch_x):
		out = []
		for i,x in enumerate(batch_x):
			valid = torch.nonzero(torch.sum(x != 0, dim=1)).squeeze()#padding are lines with all zeros	
			x = x[valid]#revome the padding
			p, t = x[:,:3], x[:,3:]
			trans = self.transform1(p)
			p = torch.mm(p, trans)
			p = torch.relu(self.fc1(p))
			p = torch.relu(self.fc2(p))
			x = torch.cat((p,t), dim=1)
			x  = torch.mm(x, self.transform2(x))
			x = torch.relu(self.fc3(x))
			x = torch.relu(self.fc4(x))
			x, _ = torch.max(x, dim=0)
			x = torch.relu(self.fc5(x))
			x = torch.relu(self.fc6(x))
			x = self.fc7(x)
			out.append(x)
		return torch.stack(out)
		

class Point_dataset(Dataset):
	def __init__(self, numpy_x, numpy_y):
		assert len(numpy_x) == len(numpy_y)
		self.x = Tensor(numpy_x)
		self.y = Tensor(numpy_y)
		
	def __getitem__(self, index):
		return self.x[index], self.y[index]
		
	def __len__(self):
		return len(self.x)










