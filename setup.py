
from setuptools import setup, find_packages

setup(
	name = "roadmap bpmatch",
	description = "",
	author = "Kyle Treleaven",
	author_email = "ktreleav@gmail.com",
	version = "0.0.0",
	packages = find_packages(),
	namespace_packages = [ 'setiptah', 'setiptah.roadbm', 'setiptah.nxopt', 'setiptah.vehrouting', ],
	install_requires = ['numpy', 'scipy', 'networkx', 'bintrees'],
)

