#Usage: cctbx.python tls_ensemble_generator_working.py pdb=$$$.pdb sampling=4 sigma=2 selection=parallel
#combined_models.pdb is the sequential TLS ensemble!
#Now, we just needt to get the probabilities
from iotbx import pdb
import sys
import math
import mmtbx.tls.tools
import tls_as_xyz_sequential
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import miller
from iotbx import scalepack

def run(args):

	vars = getargs(args)

	unit_cell, tls_segments, atom_files = extract_tls_information(vars['pdb'])

	tls_sampling = int(args[1].split('=')[1])
	resolution = float(args[2].split('=')[1])
	sigma = int(args[3].split('=')[1])
	
	tls_pdb_ensembles, sigma_info = sequential_ensemble_generator(unit_cell, tls_segments, atom_files, vars['sampling'], vars['sigma'])


	if vars['selection'] == "parallel":

		model_final = model_combining(vars['pdb'], tls_pdb_ensembles)
		get_probability('combined_models.pdb', unit_cell, sigma_info, vars['sampling'], vars['sigma'], 'parallel')



	if vars['selection'] == "anti":
		model_reverse(tls_pdb_ensembles)
		model_final = model_combining(vars['pdb'], tls_pdb_ensembles)
		get_probability('combined_models.pdb', unit_cell, sigma_info, vars['sampling'], vars['sigma'], 'anti')


def getargs(args):
	info = dict()

	for arg in args:
		data = arg.split('=')
		if len(data) == 2:
			info[data[0]]=data[1]

	info['sampling'] = int(info['sampling'])
	info['resolution'] = float(info['resolution'])
	info['sigma'] = int(info['sigma'])

	return info

def extract_tls_information(file):

	pdb_inp = pdb.input(file_name=file)
  	h = pdb_inp.construct_hierarchy()

  	unit_cell = pdb_inp.crystal_symmetry_from_cryst1()

  	pdb_obj = pdb.hierarchy.input(file_name=file)
  	remark_3_records = pdb_inp.extract_remark_iii_records(3)
  	tls_extract = mmtbx.tls.tools.tls_from_pdb_inp(
    	remark_3_records = remark_3_records, 
    	pdb_hierarchy = h)

  	tls_segments = list()
  	tls_objects = list()
  	for i, tls_object in enumerate(tls_extract.tls_params):
  		tls_segments.append(tls_object.selection_string)
  		tls_objects.append(tls_object)
  		print i, tls_object

  	count = 0

  	atom_files = list()

	for group in tls_segments:

		sel = h.atom_selection_cache().selection(string = group)
		h_new = h.select(sel)

		#Write atomic coordinates of each 
		h_new.write_pdb_file(file_name='model_%d.pdb' %count)

		atom_files.append('model_%d.pdb' %count)

		count += 1

	#Passing list of PDB files for each TLS group, unit cell parameters and TLS objects
	return unit_cell, tls_objects, atom_files


def sequential_ensemble_generator(unit_cell, tls_segments, atom_files, sample_size, sigma):

	length = len(atom_files)

	ensemble_files = list()

	sigma_info = dict()

	for i in range(0,length):
		sigma_info[i] = dict()
		#tls_as_xyz_sequential needs to return sig_x, sig_y, sig_z... for each
		sigma_info[i]['sigmas']=tls_as_xyz_sequential.run(pdb_file_name = atom_files[i], tls_object=tls_segments[i], crystal_info=unit_cell, n_models=sample_size, sigma=sigma, output_file_name='model_%d_tls.pdb' %i)
		ensemble_files.append('model_%d_tls.pdb' %i)

	return ensemble_files, sigma_info

def model_reverse(files):

	alternate = 0

	for file in files:
		if alternate%2 == 0:

			pdb_obj = pdb.input(file_name=file)
			h = pdb_obj.construct_hierarchy()

			models = []
			n_mod = 0

			for m in h.models():
				models.append(m)
				n_mod += 1


			#Define new PDB hierarchy
			output = pdb.hierarchy.ext.root()

			count = 0

			for i in range(0,n_mod):

				m_num = n_mod-count-1

				mod = models[m_num].detached_copy()
				#c = n_mod-count-1
				mod.id = models[count].id

				output.append_model(mod)

				count += 1


			output.write_pdb_file(file_name='output.pdb')


			fin = open('output.pdb', 'r')
			fout = open(file, 'w')

			lines = fin.readlines()

			for line in lines:
				data = line.split()
				if data[0] != 'BREAK':
					fout.write(line)

			fin.close()
			fout.close()
		alternate += 1

	return

def model_combining(original, models):
	data = dict()

	data['original'] = dict()
	data['original']['pdb object'] = pdb.input(file_name=original)
	data['original']['hierarchy'] = data['original']['pdb object'].construct_hierarchy()
	data['original']['models'] = []
	data['original']['chains'] = []

	for model in data['original']['hierarchy'].models():
		data['original']['models'].append(model.id)

	for chain in data['original']['hierarchy'].chains():
		data['original']['chains'].append(chain.id)


	chains = []
	for item in data['original']['chains']:
		if item not in chains:
			chains.append(item)

	#Populating the dictionary with necessary info about each TLS model
	for model in models:

		data[model] = dict()
		data[model]['pdb object'] = pdb.input(file_name=model)
		data[model]['hierarchy'] = data[model]['pdb object'].construct_hierarchy()
		data[model]['models'] = []

		count = 0
		for model_h in data[model]['hierarchy'].models():
			data[model]['models'].append(model_h)
			count += 1

		data['Number of models'] = count

	#Define the new PDB hierarchy to be written to
	output = pdb.hierarchy.ext.root()

	m_count = 0

	for i in range(0,data['Number of models']):
		new_model = pdb.hierarchy.ext.model()
		new_model.id = str(m_count)

	#Populate the new model with the chains present in the original model
		chain_list = []
		for item in chains:
			new_chain = pdb.hierarchy.ext.chain()
			new_chain.id = item
			chain_list.append(new_chain)

		for model in models:
			m_working = data[model]['models'][m_count]
			m_chains = m_working.chains()

			for ch in m_chains:
				name = ch.id

				for nm_c in chain_list:
					nm_name = nm_c.id

					if name == nm_name:
						resi = ch.residue_groups()

						for r in resi:
							re = r.detached_copy()
							nm_c.append_residue_group(re)

	#Once all the atoms across all the models are appended and accounted for, attach all the chains to the new model
		for nm_c in chain_list:
			new_model.append_chain(nm_c)

		output.append_model(new_model)
		m_count += 1



	output.write_pdb_file(file_name='output.pdb')

	fin = open('output.pdb', 'r')
	fout = open('combined_models.pdb', 'w')

	lines = fin.readlines()

	for line in lines:
		data = line.split()
		if data[0] != 'BREAK':
			fout.write(line)

	fin.close()
	fout.close()

	return 'combined_models.pdb'


def get_probability(pdb_files, unit_cell, sigma_info, tls_sampling, tls_sigma, choice):

	def s_p(sig,sig_fraction):

		probability = math.exp(-1*sig_fraction*sig_fraction/2)/(math.sqrt(2*math.pi)*sig)

		return probability

	count = 0

	#Each key in sigma_info corresponds to a specific TLS ensemble
	for key in sigma_info:
		probability = list()
		for i in range(0,tls_sampling):
			x = 2*(tls_sigma*i/(tls_sampling-1)) - tls_sigma

			prob = 1
			for s_val in sigma_info[key]['sigmas']:
				if s_val > 1.e-7:
					prob *= s_p(s_val,x) 

			probability.append(prob)

		scale_factor = 0

		for item in probability:
			scale_factor+= item

		final_probability = list()

		for item in probability:
			x = item/scale_factor
			final_probability.append(x)


		sigma_info[key]['probabilities'] = final_probability


	final_probabilities = []

	print sigma_info

	if choice == 'parallel':

		t_s = int(tls_sampling)
		for i in range(0,t_s):
			for key in sigma_info:
				probability = 1
				probability *= sigma_info[key]['probabilities'][i]
			final_probabilities.append(probability)


	if choice == 'anti':

		t_s = int(tls_sampling)
		for i in range(0,t_s):
			for key in sigma_info:
				probability = 1

				if key%2 != 0:
					probability *= sigma_info[key]['probabilities'][tls_sampling-i-1]

				else:
					probability *= sigma_info[key]['probabilities'][i]
			final_probabilities.append(probability)

	scale_factor = 0

	for item in final_probabilities:
		scale_factor += item


	actual_probs = []
	sum = 0
	for item in final_probabilities:
		x = item/scale_factor

		actual_probs.append(x)
        print '***HERE ARE THE MODEL PROBABILITIES!!!***'
	print actual_probs


def map_combining(files):
	output = dict()

	count = 0

	for file in files:
		fin = open(file, 'r')
		lines = fin.readlines()

		for line in lines:

		
			data = line.split()

			h = int(data[0])
			k = int(data[1])
			l = int(data[2])
			signal = data[3]

			
			signal_new = float(signal)

			if h not in output:
				output[h] = dict()
			if k not in output[h]:
				output[h][k] = dict()
			if l not in output[h][k]:
				output[h][k][l] = 0
		
			output[h][k][l] += signal_new/len(files)

		fin.close()
		count += 1

	fout = open('final_result.hkl', 'w')
	for key_h in output:
		for key_k in output[key_h]:
			for key_l in output[key_h][key_k]:
				print >>fout, "%4d %4d %4d %4d" %(key_h, key_k, key_l, output[key_h][key_k][key_l])

	fout.close()

	print 'Everything ran successfully!'

	return 'final_result.hkl'

if __name__ == "__main__":
	import sys
	run(sys.argv[1:])
