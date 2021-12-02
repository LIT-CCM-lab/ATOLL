import sys
import os
from copy import deepcopy

import sequence
import frame
import procedure
import plot
import initialize as init

__authors__ = ('Célien Jacquemard', 'Guillaume Bret')
__version__ = (0, 1)
__creation_date__ = '2020.07.16'
__last_update__ = '2020.07.16'
__email__ = ('jacquemard@unistra.fr' , 'gbret@unistra.fr')
__status__ = 'InDev'

prog_dir = os.path.realpath(os.path.dirname(__file__))
logger = init.get_logger(__file__)

class Atoll:
	def __init__(self):
		self.reference = None
		self.structures = {}

	def load_from_cli(self, args):
		# Load sequence file
		if args.sequence_filepath:
			self.sequence_loader = sequence.SequenceLoader(args.sequence_filepath)
		else:
			self.sequence_loader = sequence.BasicSequenceLoader()

		# Load structure informations
		if args.info_filepath:
			self.info = frame.Info(args.info_filepath)
		else:
			self.info = frame.BasicInfo()

		# Domain definition
		self.domains = {}
		self.resnum = args.resnum
		logger.info(f'Numbering scheme: {args.resnum}')

		# Residues use for reference centering
		if args.center:
			self.domains['center'] = frame.Selector.split(args.center)

		# Residues use for reference principal axes
		if args.pmi:
			self.domains['pmi'] = frame.Selector.split(args.pmi)

		if args.pmi_x and args.pmi_y and args.pmi_z:
			self.pmi_axes = args.pmi_x, args.pmi_y, args.pmi_z
		else:
			self.pmi_axes = (2, 0, 1)

		# Residues use for structure alignment
		if args.resalign:
			self.domains['alignment'] = frame.Selector.split(args.resalign)

		# Residues use for analysis
		if args.reshelix:
			self.domains['phelix'] = frame.Selector.split(args.reshelix)
		else:
			raise NotImplementedError

		# Load reference and query structures
		self.reference = self.load_reference(args.reference_filepath)
		self.structures = self.load_structures(self.info, args.structure_dirpath)

		if not self.structures:
			logger.critical('No input files were loaded!')
			return 1

		# Load protein informations
		self.set_info(self.reference)
		for structure in self.structures:
			self.set_info(structure)

		# Define amino acid sequences for reference and structures
		self.set_sequence(self.reference)
		for structure in self.structures:
			self.set_sequence(structure)

		# Set domains to reference and queries
		self.set_domain(self.reference)
		for structure in self.structures:
			self.set_domain(structure)

		# Align reference sequence onto structure sequence
		sequence_aligner = sequence.SequenceAligner()
		sequence_aligner.align_reference(self.reference)
		for structure in self.structures:
			logger.info(f'Sequence alignment of "{structure.label}"')
			sequence_aligner.align_reference(structure)

		self.output = frame.Output(args.output_dirpath, args.overwrite)
		if args.center is not None:
			center_atoms = self.reference.get_residues_by_group('center').atoms.select_atoms('name CA')
			self.reference.center(center_atoms)

		if args.pmi is not None:
			pmi_atoms = self.reference.get_residues_by_group('pmi').atoms.select_atoms('name CA')
			self.reference.align_pmi(pmi_atoms, self.pmi_axes)

		if args.transformation_matrix is not None:
			try:
				tm = [float(x) for x in args.transformation_matrix.strip('"').strip("'").split('|')]
			except ValueError as error:
				logger.critical(f'Your transformation matrix is not properly formatted {args.transformation_matrix}')
				return 1
			logger.info('Transform with matrix:')
			logger.info(tm)
			self.reference.transform(tm)

		self.reference.set_referential()

		self.analyzer = procedure.Analyzer(self.structures, self.reference, self.output.dirpath)
		if args.resalign is not None:
			self.analyzer.add_analysis('align')
		# self.analyzer.add_analysis('trajwrite')
		self.analyzer.add_analysis('phelix')
		self.analyzer.run()

		pm = plot.ProjectionMap(
			self.analyzer.analysis_procedures[-1],
			output_filepath=os.path.join(self.output.plot_dirpath, 'projection_map.'+args.image_format),
			merging_type=args.merging_type,
			colors=self.info.get_dict_values('color')
		)
		logger.info('DONE')

		return 0

	@classmethod
	def load_reference(cls, path, *args, **kwargs):
		logger.info(f'Loading reference {os.path.basename(path)}...')
		reference = frame.Reference(path, *args, **kwargs)
		return reference

	@classmethod
	def load_structures(cls, info, base_dirpath=None):
		structures = []
		for entry in info.entries.values():
			if base_dirpath:
				path = os.path.join(base_dirpath, entry.path)
			else:
				path = entry.path
			label = entry.id
			fileprefix, fileext = os.path.splitext(os.path.basename(path))

			if os.path.isfile(path):
				structure = frame.Static(label, path)
			else:
				structure = frame.Multiple(label, path)

			structures.append(structure)

		return structures

	def set_info(self, frame):
		if frame.label in self.info:
			frame.load_info(self.info[frame.label])

	def set_sequence(self, frame):
		if frame.label in self.sequence_loader:
			logger.info('by label')
			frame_sequence = self.sequence_loader[frame.label]
		elif frame.protein_name in self.sequence_loader:
			logger.info(f'by protein name {frame.protein_name}')
			frame_sequence = self.sequence_loader[frame.protein_name]
		elif self.sequence_loader.reference is not None:
			logger.info('by reference')
			frame_sequence = self.sequence_loader.reference
		elif self.reference:
			logger.info('by reference structure')
			frame_sequence = sequence.Sequence(self.reference.protein, fill_empty=False, position_shift=False)
		else:
			logger.info('by entry structure')
			frame_sequence = sequence.Sequence(frame.protein, fill_empty=False, position_shift=False)

		frame.sequence = deepcopy(frame_sequence)

	def set_domain(self, frame):
		logger.info(f'Setting domains for {frame.label}...')
		for domain_name, intervals in self.domains.items():
			frame.domains.set_group_intervals(domain_name, intervals, self.resnum, prefix='TM')


class RefNum(Atoll):
	def __init__(self):
		self.reference = None
		self.sequence = None

	def load_from_cli(self, args):

		if args.sequence_filepath:
			self.reference = self.load_reference(args.reference_filepath, fill_empty=False)
			self.sequence_loader = sequence.SequenceLoader(args.sequence_filepath)
			self.set_sequence(self.reference)

			sequence_aligner = sequence.SequenceAligner()
			new_sequence, new_protein = sequence_aligner.reference_renumber(self.reference)
		else:
			self.reference = self.load_reference(args.reference_filepath, fill_empty=False, position_shift=False)
			self.reference.structure_sequence.left_strip()

			new_protein = self.reference.protein
			new_sequence = self.reference.structure_sequence

		import csv
		with open(args.output_filepath, 'w') as f:
			table = csv.writer(f, delimiter='\t')
			for aa in new_sequence:
				table.writerow([aa.tl, aa.ol, aa.index, aa.position, aa.resid])

		if args.output_reference_filepath:
			new_protein.write(args.output_reference_filepath)


class CheckFile(Atoll):
	def __init__(self):
		pass

	@staticmethod
	def format_error(error):
		return ' '.join(str(error).split('\n'))

	def load_from_cli(self, args):
		import csv
		reports = []

		row = ['Sequence', args.sequence_filepath]
		if args.sequence_filepath:
			try:
				sequence_loader = sequence.SequenceLoader(args.sequence_filepath)
				row.append('Ok')
				row.append(None)
			except Exception as error:
				row.append('Failed')
				row.append(CheckFile.format_error(error))
		else:
			row.append('Missing')
			row.append(None)
		reports.append(row)

		row = ['Annotation', args.info_filepath]
		info = None
		if args.info_filepath:
			try:
				info = frame.Info(args.info_filepath)
				row.append('Ok')
				row.append(None)
			except BaseException as error:
				row.append('Failed')
				row.append(CheckFile.format_error(error))
		else:
			row.append('Missing')
			row.append(None)
		reports.append(row)

		row = ['Reference', args.reference_filepath]
		if args.reference_filepath:
			try:
				reference = self.load_reference(args.reference_filepath)
				row.append('Ok')
				row.append(None)
			except Exception as error:
				row.append('Failed')
				row.append(CheckFile.format_error(error))
		else:
			row.append('Missing')
			row.append(None)
		reports.append(row)

		row = ['Structure', args.info_filepath]
		if args.structure_filepathes:
			try:
				if info is None:
					for label, path in enumerate(args.structure_filepathes):
						if os.path.isfile(path):
							structure = frame.Static(label, path)
						else:
							structure = frame.Multiple(label, path)
				else:
					structures = info.load_structures()
				row.append('Ok')
				row.append(None)
			except Exception as error:
				row.append('Failed')
				row.append(CheckFile.format_error(error))
		else:
			row.append('Missing')
			row.append(None)
		reports.append(row)

		with open(args.report_filepath, 'w') as f:
			table = csv.writer(f, delimiter='\t')
			for row in reports:
				table.writerow(row)


class Sanitize(Atoll):
	def __init__(self):
		pass

	def load_from_cli(self, args):
		structure = frame.Static('multimeric', args.structure_filepath)
		backbone = structure.protein.select_atoms('backbone')
		for i, residue in enumerate(backbone.residues, 1):
			residue.resid = i
			residue.segment.segid = 'A'

		backbone.write(args.output_filepath)



def main(args):
	atoll = Atoll()
	return atoll.load_from_cli(args)


def refnum(args):
	refnum = RefNum()
	refnum.load_from_cli(args)

def checkfile(args):
	checkfile = CheckFile()
	checkfile.load_from_cli(args)

def sanitize(args):
	sanitize = Sanitize()
	sanitize.load_from_cli(args)

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(
		prog='atoll',
		description='The ATOLL program analyzes transmembrane protein helix bundle.',
		epilog='For more details, see README.md.'
	)

	subparsers = parser.add_subparsers()
	atoll_parser = subparsers.add_parser('atoll', help='')

	input_group = atoll_parser.add_argument_group('input')
	input_group.add_argument('--reference', '-ref', required=True, dest='reference_filepath',
		help='Reference structure file.', metavar='FILEPATH')
	input_group.add_argument('--sequence', '-seq', dest='sequence_filepath',
		help='Sequence file of studied proteins.', metavar='FILEPATH')
	input_group.add_argument('--info', '-inf', dest='info_filepath', required=True,
		help='Information file for each entry.', metavar='FILEPATH')
	input_group.add_argument('--structure', '-s', dest='structure_dirpath',
		help='Information file for each entry.', metavar='DIRPATH')

	selection_group = atoll_parser.add_argument_group('selection')
	selection_group.add_argument('--resnum', '-rn', choices=['position', 'resid'], default='position',
		help='Residue numbering to apply.')
	selection_group.add_argument('--resalign', '-ra',
		help='Selection of residues involved in structure alignment.')
	selection_group.add_argument('--reshelix', '-rh', required=True,
		help='Selection of residues to analyze that correspond to transmembrane helice ends.')
	selection_group.add_argument('--center',
		help='Selection of residues that will be centered on the origin.')
	selection_group.add_argument('--pmi',
		help='Selection of residues whose principal axes of inertia are aligned with the axes of the reference frame.')
	selection_group.add_argument('--pmi_x', metavar='X', type=int,
		help='')
	selection_group.add_argument('--pmi_y', metavar='Y', type=int,
		help='')
	selection_group.add_argument('--pmi_z', metavar='Z', type=int,
		help='')
	selection_group.add_argument('--tmatrix', '-tm', dest='transformation_matrix',
		help='The transformation matrix applied to the reference structure')

	output_group = atoll_parser.add_argument_group('output')
	output_group.add_argument('--output', '-out', required=True, dest='output_dirpath',
		help='Output directory where results will be stored.', metavar='DIRPATH')
	output_group.add_argument('--overwrite', action='store_true',
		help='Overwrite data if output directory is already exist.')

	phelix_group = atoll_parser.add_argument_group('atoll')
	phelix_group.add_argument('--merge', default=None, choices=[None, 'class'], dest='merging_type',
		help='Merge contour with same class name.', metavar='MERGING_TYPE')
	phelix_group.add_argument('--image_format', default='png', choices=['png', 'svg'], dest='image_format',
		help='Projection image format.', metavar='FILE_FORMAT')

	misc_group = atoll_parser.add_argument_group('miscellaneous')
	# Optional argument that count the number of flag occurence (-v or -vv).
	misc_group.add_argument('-v', '--verbosity', action='count', default=0,
		help='Increase output verbosity.')

	atoll_parser.set_defaults(func=main)

	refnum_parser = subparsers.add_parser('refnum', help='Renumber the residue in a structure file based on sequence alignment file.')
	refnum_parser.add_argument('--reference', '-ref', required=True, dest='reference_filepath',
			help='Reference structure file.', metavar='FILEPATH')
	refnum_parser.add_argument('--sequence', '-seq', dest='sequence_filepath',
		help='Sequence file of studied proteins.', metavar='FILEPATH')
	refnum_parser.add_argument('--output', '-o', dest='output_filepath', required=True,
		help='The output file containing residue indices and positions.', metavar='FILEPATH')
	refnum_parser.add_argument('--output_reference', '-or', dest='output_reference_filepath',
		help='The output reference with protein only.', metavar='FILEPATH')

	refnum_parser.set_defaults(func=refnum)

	checkfile_parser = subparsers.add_parser('checkfile', help='Module to check if input files are readable by ATOLL.')
	checkfile_parser.add_argument('--reference', '-ref', dest='reference_filepath',
		help='Reference structure file.', metavar='FILEPATH')
	checkfile_parser.add_argument('--sequence', '-seq', dest='sequence_filepath',
		help='Sequence file of studied proteins.', metavar='FILEPATH')
	checkfile_parser.add_argument('--info', '-inf', dest='info_filepath',
		help='Information file for each entry.', metavar='FILEPATH')
	checkfile_parser.add_argument('--structure', '-struc', dest='structure_filepathes', nargs='*',
		help='Structure entries.', metavar='FILEPATH or DIRPATH')
	checkfile_parser.add_argument('--output', '-o', dest='report_filepath',
		help='The output report file.', metavar='FILEPATH')

	checkfile_parser.set_defaults(func=checkfile)

	mergechains_parser = subparsers.add_parser('sanitize', help='Clean a input structure to suit ATOLL procedure.')
	mergechains_parser.add_argument('--structure', '-s', dest='structure_filepath',
		help='Structure file describing only ONE conformer.', metavar='FILEPATH')
	mergechains_parser.add_argument('--output', '-o', dest='output_filepath',
		help='', metavar='FILEPATH')

	mergechains_parser.set_defaults(func=sanitize)

	args = parser.parse_args()

	# Run
	status = args.func(args)
	sys.exit(status)
