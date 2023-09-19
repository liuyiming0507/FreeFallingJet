# FreeFallingJet

The mesh files are too big to upload, so please find the whole geofolder in google drive
https://drive.google.com/drive/folders/129dJLroMg20GUEs7fWGnHnS-RLbWxP1z?usp=sharing

## FEniCS
**These codes base on the *FEniCS* not *FEniCS-X*.**
### Fast way of building FEniCS on Docker (WINDOWS/MACOS/UBUNTU)
1. Download Docker Desktop: https://www.docker.com/products/docker-desktop/
2. Run Docker Desktop and open the terminal to run `docker run hello-world` for checking that Docker is working.
3. Access the folder where the local code is located (cd command or directly open the folder and open the terminal through the menu of the right mouse button in the folder).
4. Run `docker run -ti -v $(pwd):/home/fenics/shared quay.io/fenicsproject/stable` in MacOS and Ubuntu, <br>
   or run  `docker run -ti -v ${pwd}:/home/fenics/shared quay.io/fenicsproject/stable` in Windows. <br>
   The difference is `$(pwd)` or `${pwd}`.
5. Run `ls` to list the code files under this shared folder and you can find the code files in your local folder here.
6. Run `mpirun -n 8 python3 yourFEniCScodeName.py` or whatever commend you like.

--------------***Step 4 is critical***--------------

**Reference here https://fenics.readthedocs.io/projects/containers/en/latest/introduction.html#installing-docker**
