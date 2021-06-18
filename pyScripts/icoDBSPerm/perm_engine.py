import re
import os, shutil
import numpy as np
import logging
from logger import LoggerService
from of_runner import OFRunner
from am_domain_reader import AMDomainReader

class PermEngine(object):
    def __init__(self,domain_reader,template_config) -> None:
        super().__init__()
        self._domain_reader=domain_reader
        self._geometry_template_dir=template_config["geometry"]
        self._perm_case_template_dir=template_config["perm_case"]

        numeric_const_pattern = r"""
            [-+]? # optional sign
            (?:
                (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
                |
                (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
            )
            # followed by optional exponent part if desired
            (?: [Ee] [+-]? \d+ ) ?
            """
        self._numeric_detecter = re.compile(numeric_const_pattern, re.VERBOSE)
    
    def run(self,case):
        case_dir=case["directory"]
        structure_file=case["structure"]
        nx=case["nx"]
        ny=case["ny"]
        nz=case["nz"]
        pixel_resolution=case["resolution"]
        subdomain_num=case["subdomain_num"]
        np=case["processor"]
        
        logger_service = LoggerService(case_dir,logging.DEBUG)
        case_logger=logger_service.get_logger()

        case_logger.debug("about to read the domain file")
        domain=self._domain_reader.read(structure_file,nx,ny,nz)
        case_logger.info(f"domain shape: {domain.shape}")

        case_logger.debug("about to chunk the domain")
        nz_chunks=self._chunk(range(nz), subdomain_num)
        subdomains=[]
        for chunk in nz_chunks:
            subdomain=domain[:,:,chunk]
            subdomains.append(subdomain)

        case_logger.debug("about to create each subdomain case")
        subdomain_case_dirs=[]
        for i in range(subdomain_num):
            subdomain_dir=os.path.abspath(os.path.join(case_dir,f"subdomain{str(i)}"))
            subdomain_case_dirs.append(subdomain_dir)
            self._copytree(self._perm_case_template_dir, subdomain_dir)
            self._create_geometry_files(subdomains[i],self._geometry_template_dir,subdomain_dir,pixel_resolution)
            self._modify_computation_processor_num(subdomain_dir,np)
            case_logger.debug(f"succeed creating subdomain case: {subdomain_dir}")
        
        case_logger.info(f"Succeed creating {subdomain_num} subdomain cases")

        case_logger.debug("about to compute the subdomain permeabilities")
        of_runner=OFRunner(logger_service)
        subdomain_perms=[]
        try:
            for subdomain,subdomain_case_dir in zip(subdomains,subdomain_case_dirs):
                case_logger.info(f"running submain case: {subdomain_case_dir}")
                subdomain_name=os.path.basename(subdomain_case_dir)
                porosity=self._report_porosity(subdomain)
                case_logger.info(f"{subdomain_name} porosity: {porosity}")
                perm=of_runner.run(subdomain_case_dir)
                case_logger.info(f"{subdomain_name} permeability: {perm}")
                subdomain_perms.append(perm)
        except Exception as e:
            err_msg=f"Exception happended when running the openFoam runner, with error message: \n {e} \n" + \
                    "Check run.log for more details"
            case_logger.error(err_msg)
            return -1,False,err_msg

        
        case_logger.debug("about to calculate the whole domain permeability")
        subdomain_nzs=[]
        for subdomain_case_dir in subdomain_case_dirs:
            blockMesh_file=os.path.join(subdomain_case_dir,"system/blockMeshDict")
            with open(blockMesh_file,"r") as fp:
                lines=fp.readlines()
            for line in lines:
                if line.startswith("Nz"):
                    subdomain_nz=int(self._numeric_detecter.findall(line)[0])
                    subdomain_nzs.append(subdomain_nz)
                    break
        
        perm=0
        for subdomain_nz,subdomain_perm in zip(subdomain_nzs,subdomain_perms):
            perm+=subdomain_nz/subdomain_perm
        perm=nz/perm

        case_logger.info(f"Domain Permeability:{perm}")
        return perm,True,""


    def _chunk(self,seq, num):
        avg = len(seq) / float(num)
        out = []
        last = 0.0

        while last < len(seq):
            out.append(seq[int(last):int(last + avg)])
            last += avg

        return out 

    def _report_porosity(self,domain):
        label_pore=0
        # label_rock=1
        porosity=np.sum(domain==label_pore)/domain.size
        return porosity

    def _copytree(self,src, dst, symlinks=False, ignore=None):
        if os.path.exists(dst):
            shutil.rmtree(dst)

        for item in os.listdir(src):
            s = os.path.join(src, item)
            d = os.path.join(dst, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, symlinks, ignore)
            else:
                shutil.copy2(s, d)

    def _create_geometry_files(self,subdomain,geometry_template_dir,subdomain_dir,pixel_resolution):
        nx=subdomain.shape[0]
        ny=subdomain.shape[1]
        nz=subdomain.shape[2]
        solid_internal_field=list()
        solid_internal_field.append(f"{str(nx*ny*nz)}\n")
        solid_internal_field.append("(\n")
        for k in np.arange(0,nz):
            for j in np.arange(0,ny):
                for i in np.arange(0,nx):
                    solid_internal_field.append(f"{str(subdomain[i,j,k])}\n")
        solid_internal_field.append(")\n")
        solid_internal_field.append(";\n")

        with open(f"{geometry_template_dir}/solid_template","r") as fp:
            solid_template=fp.readlines()
        internal_field_line_index=0
        for index,line in enumerate(solid_template):
            if line.startswith("internalField"):
                internal_field_line_index=index
                for i in range(len(solid_internal_field)):
                    solid_template.insert(internal_field_line_index+1+i,solid_internal_field[i])
                break

        if not os.path.exists(f"{subdomain_dir}/0"):
            os.makedirs(f"{subdomain_dir}/0")
        
        if not os.path.exists(f"{subdomain_dir}/system"):
            os.makedirs(f"{subdomain_dir}/system")

        with open(f"{subdomain_dir}/0/solid","w") as fp:
            fp.writelines(solid_template)

        with open(f"{geometry_template_dir}/blockMeshDict_template","r") as fp:
            mesh_template=fp.readlines()
        for index,line in enumerate(mesh_template):
            if line.startswith(r"//insert here"):
                mesh_template.insert(index+1,f"convertToMeters {pixel_resolution};\n")
                mesh_template.insert(index+2,f"Nx {nx};\n")
                mesh_template.insert(index+3,f"Ny {ny};\n")
                mesh_template.insert(index+4,f"Nz {nz};\n")
                break

        with open(f"{subdomain_dir}/system/blockMeshDict","w") as fp:
            fp.writelines(mesh_template)

    def _modify_computation_processor_num(self,subdomain_dir,np):
        file_dir=f"{subdomain_dir}/system/decomposeParDict"
        with open(file_dir,"r") as fp:
            decomposeParDict=fp.readlines()
        for index,line in enumerate(decomposeParDict):
            if line.startswith("numberOfSubdomains"):
                decomposeParDict[index]=f"numberOfSubdomains {np};\n"
                break
        with open(file_dir,"w") as fp:
            fp.writelines(decomposeParDict)


if __name__ == "__main__":
    case={}
    case["directory"]="/home/anoldfriend/OpenFOAM/anoldfriend-7/applications/solvers/icoDBSPermFoam/runs/splitcases60x50x100"
    case["structure"]="/home/anoldfriend/OpenFOAM/anoldfriend-7/applications/solvers/icoDBSPermFoam/runs/splitcases60x50x100/structures/DATA_60x50x100.am"
    case["nx"]=60
    case["ny"]=50
    case["nz"]=100
    case["resolution"]=25e-6
    case["subdomain_num"]=2
    case["processor"]=8

    template_config={}
    template_config["geometry"]="/home/anoldfriend/OpenFOAM/anoldfriend-7/applications/solvers/icoDBSPermFoam/runs/splitcases60x50x100/geometry_template"
    template_config["perm_case"]="/home/anoldfriend/OpenFOAM/anoldfriend-7/applications/solvers/icoDBSPermFoam/runs/splitcases60x50x100/perm_case_template"

    
    reader=AMDomainReader()
    perm_engine=PermEngine(reader,template_config)
    perm_engine.run(case)
