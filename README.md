# Simulation Best Practices for Lipid Membranes

Best Practices document to be submitted to the Living Journal of Computational Molecular Science

Contributors: David Smith (davesmith4398), Alan Grossfield (agrossfield), Ed Lyman, Michael Shirts (mrshirts)
People to (possibly) invite: Scott Feller, Richard Pastor

The GitHub document is:

https://github.com/davesmith4398/best_practice_membranes


# Introduction and Scope

Goal: establish reliable and robust standardization of settings and setup practices for practical molecular dynamics simulations of pure lipid bilayer membranes.

Key terms?

Background

Lipid bilayer membranes (planar and vesicular) have applications in soft matter physics, pharmacology, consumer products, etc., and are first approximants to biological membranes. Lipid bilayers consist of two molecularly-thick layers, or leaflets, in aqueous solvent, where the lipids are oriented with their polar head groups outward and their nonpolar tail groups inward. The dominant driving forces in their formation and stability include hydrophobic and dispersion interactions, wherein the lipid tail groups maximize their contacts with each other and minimize their contacts with water, decreasing area per lipid; and head group electrostatic repulsion and tail group conformational entropy, which work to increase the area per lipid. “Fluid” (liquid crystalline) lipid membranes are normally modeled as liquid-like laterally (no in-plane shear modulus), and solid-like transversally (out of plane). In lipid membrane simulations, the canonical system is DPPC (dipalmitoylphosphatidylcholine) in water. This is typically the system for which new force fields are first tested. DPPC is sometimes not preferred in experiments, however, due to its high melting point (from the gel to liquid-crystalline state). In any event, the main goal for a lipid membrane model study should be to correctly capture the correct structural, mechanical, thermodynamic, and/or dynamic properties (whatever is relevant) at the relevant length and timescales and the correct equilibrium (thermodynamic, temperature, pressure, etc.) and/or nonequilibrium (thermal/mechanical/chemical/other gradients) conditions. Transferability is ideal, but not necessarily a priority.

Some good papers/textbooks for statistical mechanical/thermodynamic background on membranes
 - Safran, Samuel A. “Statistical Thermodynamics of Surfaces, Interfaces, and Membranes.” 2003: Westview Press.
 - Nelson, D.; Piran, T.; and Weinberg, T. “Statistical Mechanics of Membranes and Surfaces.” 2004: World Scientific Publishing Company.
 - Boal, David. “Mechanics of the Cell.” 2012: Cambridge University Press, New York.
 
Good papers/textbooks for computational/simulation guidance on membranes
 - Sundararajan, V. “Computational Modeling of Membrane Bilayers, Volume 60 (Current Topics in Membranes).” 2008: Academic Press.
 - Tieleman, Marrink, and Berendsen. Biochimica et Biophysica Acta, 1997.


# The Checklist
# I. Force field/model selection

As with other systems, model selection for lipid membranes is crucial. Since the efficiency of planar membrane simulations scales geometrically with the square of the in-plane box length and in equilibration with the fourth power of the box length, the total scaling can be up to the sixth power of the box length. Again, the main goal is accuracy in properties at the relevant spatiotemporal resolution and equilibrium/nonequilibrium resolution. As for model resolution, lipid force fields range from all-atom to united atom to coarse-grained (CG) and from explicit to implicit solvent. Atomistic models can have hundreds of atomic sites per lipid molecule, while models like the Martini CG force field can have around ten pseudoatom sites and dissipative particle dynamics (DPD) simulations can have just two or three. Because the aqueous solvent can contribute up to 90 percent of the force evaluations, implicit solvent (IS) simulations can also be very advantageous. Lower-resolution models can extend accessible length and timescales from about tens of nanometers and hundreds of nanoseconds to hundreds of nanometers and tens of microseconds, but are sometimes deficient in capturing solvent-mediated effects, entropic driving forces, and hydrodynamics. For lipid membranes, there is also an extensive sub-community that uses continuum mechanical theory and field-theoretic simulations that are sometimes in fact the preferred approach at larger length scales due to their efficiency. For lipid membrane studies, it is often crucial to compare with these techniques, and evaluate whether or not a continuum approach would be better.


***For the sections below, good to first consult the Interfacial Systems practices, which should be consistent with this one (with major differences and nuances always noted).***


# II. Preparing initial configurations

Out of convention, the membrane in-plane directions are often defined/set to be the x and y directions, while the out-of-plane direction is defined as z. All subsequent items assume this directionality, but the choice of direction is otherwise arbitrary.

Setup: Unilamellar lipid bilayers exist for a variety of concentrations. Certain force fields at a certain resolution normally prescribe some number of water molecules for every lipid. For a planar bilayer, given a desired number of lipids or membrane size in the xy-plane (which can be estimated from one another if the area per lipid is roughly known), the concentration/amount of solvent can be independently varied by changing the box z-dimension. For a vesicle, the total membrane area must be smaller than the smallest plane in periodic box; otherwise, the lipids will form a planar bilayer to minimize the free energy. There are actually two main methods for setting up planar bilayers and vesicles: (1) templating and (2) self-assembly.

 - TEMPLATING: lipid amphiphiles are pre-assembled in a planar or vesicular bilayer. Due to core overlaps for both and leaflet number asymmetry for vesicles (especially at smaller radii), existing packages and routines are probably preferred. CHARMM-GUI (http://www.charmm-gui.org/) is an excellent resource, and the most common for setting up membranes in a variety of configurations and for a variety of models/force fields. Per the general interfacial system recommendation of combining systems one at a time, the membrane can then be solvated. Due to the fluctuating nature and molecular scale roughness of lipid membranes, a trial-and-error solvation routine may be preferred over appending solvent “slabs” to each side of the membrane (although judgment can be used to determine which one will be more efficient).

    - General procedure:
        - Generate starting structure
            - CHARMM-GUI
                - Most commonly used package
                - Packs lipids from a library, then relaxes clashes
                - MUST CHECK FOR DISASTER STRUCTURES, e.g. chains going through rings, flipped chiralities of lipid backbone
            - OptimalMembraneGenerator from LOOS
                - Much less commonly used (newer, no new community)
                - NAMD only
                - Packs lipids from a library, but uses a scaling procedure like insane so no nutty structures
            - For multicomponent lipid membranes, there are some tools to help build membranes given the desired lipid types and proportions (e.g. insane.py http://www.cgmartini.nl/index.php/downloads/tools/239-insane for the Martini CG force field).
                - Insane tries to do most of your thinking for you, sacrifices “control” for correctness (won’t let you do something stupid
        - Minimize (no dynamics) to remove bad site contacts
        - Anneal from 0 K to target temperature in relevant ensemble (typically NVT or NPT, depending on whether or not you want to the box dimensions to relax/fluctuate) to slowly heat system to temperature of interest. For soft enough interactions, you may be able to skip this step.
    - Strengths
        - Directly construct a sane-looking bilayer
        - Efficient
    - Weaknesses
        - Doesn’t capture preferential segregation of multicomponent bilayers, which can be slow to emerge (so you may have a very long equilibration time that’s not obvious)
        - You need to know what you want (say what’s in what leaflet, etc), and what you want might be wrong

 - SELF-ASSEMBLY: the desired number of lipids is solvated in a water box, then allowed to spontaneously assemble into a vesicle or planar bilayer through a simulation in the NVT or NPT ensemble. The proof-of-concept approach is more scientifically satisfying, but not necessarily reasonable for higher resolution models due to the long timescale of assembly. Because the system will always proceed to minimize its total free energy, the initial box dimensions is crucial to the outcome of the self-assembly approach. For a planar bilayer, the in-plane box dimensions must be initialized near their intended end state (which can be predicted with the number of lipids and area per lipid), and for a vesicle, all box areas must be significantly larger than the corresponding planar bilayer with the same number of lipids (otherwise, the lipids will form a planar bilayer and not a vesicle). For one example of this second technique, you can visit the Martini website (http://www.cgmartini.nl/index.php/tutorials-general-introduction/bilayers).

    - General procedure:
        - Generate starting structure (lipids in solution)
        - Minimize (no dynamics) to remove bad contacts
        - Anneal from 0 K to target temperature in relevant ensemble, possibly with position restraints on all sites. For soft enough interactions, you may be able to skip this step.
        - Self-assembly in NVT or NPT. The time for this to complete can vary widely with initial configuration. Once a vesicle or planar bilayer is formed, the bilayer or vesicle can be reoriented if desired (e.g. such that x and y are the in-plane directions).
        - Equilibrate in target ensemble; equilibration times vary with systems size and phenomenon of interest
    - Strengths
        - “Natural”: lets things segment the way they want to
            - Great if you know the overall composition but not the distribution
        - Easy (at least with CG calculations) -- scatter and run
    - Weaknesses
        - Less reproducible -- bilayers won’t be symmetric, leaflets won’t necessarily have same composition.  Can be problematic with small systems
        - Relatively expensive


# III. Equilibrating the system and selecting MD settings

 - Target temperature: determined also by desired membrane phase; the main transition temperature of interest is the gel-to-liquid phase transition temperature, above which the membrane exists in a disordered liquid crystalline L$_\alpha$ state and below which the membrane exists in an ordered gel L$_\beta$ state. Some models may capture intermediate tilted gel L$_\beta’$, ripple, and interdigitated phases whose relevance depends on the experiments you are trying to model. Most simulations approximate a cellular membrane as a fluid lipid bilayer, and build heterogeneity in later.
 - Thermostat (see general interfacial system recommendations as well)
    - For CG and IS models, a stochastic dynamics integrator/thermostat (e.g. Langevin or Brownian dynamics) may be prescribed. The Langevin thermostat is not guaranteed to accurately capture long-range hydrodynamics, but efforts have been made to improve hydrodynamics for IS CG models (see: Lyman, Atzberger).
 - Target ensemble: the main ones of relevance for lipid membranes are NVT and NP$_z$P$_xy$T (or NP$_z$tT, the multiphase ensemble, where t is the frame tension). The latter involves pressure control in a semiisotropic scheme, which controls the xy and z pressures independently. Since soft matter and biological membranes often operate at negligible tension (as their conjugate variable, the area per lipid, is unconstrained and therefore used to minimize the free energy), tensionless membranes are currently the most common. NVT is appropriate for closed vesicles, due to their isotropic nature. The multiphase ensemble (semiisotropic pressure coupling) is preferred for planar bilayers. For membranes, the definition of tension is a precarious one that might not be trivial to a newcomer. While there is an important distinction between the frame tension, conjugate to the frame area (the box x and y dimensions), and the Laplace tension , conjugate to the fluctuating membrane contour area, it has been clearly shown through thermodynamic arguments that these tensions and areas are directly related, and therefore not independent (Diamant. Phys. Rev. E, 2011.). The Laplace tension is defined to a first approximation as $\gamma$ = L$_z$/2(P$_z$-(P$_x$+P$_y$)/2), and both tensions reduce to zero when the component pressures are set to be equal.
    - However, phase coexistence studies should employ NVT, as the total area for the multiphase (e.g. fluid and gel) membrane must be intermediate to the total areas of the pure fluid and pure gel membranes (with proportions of fluid and gel phase lipids determinable, e.g. via the lever rule).
    - Additional (if relevant to simulation ensemble):
        - Compressibilities: inverse to some other interfacial simulations (e.g. SAMs on a gold surface in water), where the in-plane compressibilities are set to zero to preserve hydrocarbon area per molecule, tilt, and density, membrane simulations are usually set to be compressible in the xy plane, and sometimes even incompressible in z (especially for IS CG models).

***Production run section for this?? More of a sampling issue than an equilibration one.***
 - Time: depending on the phenomenon of interest, you will ultimately want to incorporate into your production run an adequate amount of sampling. In general, this means several autocorrelation times for the relevant degrees of freedom. The time of large-wavelength membrane undulations, for example, will scale as ~L$^4$, where L is the x or y (in-plane) dimension. Undulations, interleaflet reorganization (e.g. lipid flip flop), and other collective order fluctuations, by definition, involve the coordination and motion of several lipids, and therefore will be the longest timescale fluctuations in the system (Vermeer et al. Eur. Biophys. J., 2007.). However, if your study concerns more localized or molecular degrees of freedom (e.g. rotation about chemical bonds, trans-gauche isomerization, lipid axial diffusion, etc.) at constant membrane macroscopic shape and lipid leaflet number, then smaller sampling times are permissible.


# IV. Properties to check (are you simulating what you want?)

“Fluid” (liquid crystalline) lipid membranes are normally modeled as liquid-like laterally (no in-plane shear modulus), and solid-like transversally. Because of this, important properties include in-plane elasticity/dynamics and out-of-plane elasticity.

 - Lateral density profile/membrane thickness: Various regional models have been proposed to characterize the membrane based on its component and overall densities.
    - Compute X-ray scattering (electron density)
    - Neutron scattering
    - Must compute scattering profiles and compare to experiment, NOT compare density profile to the “experimental ones”
    - Electrostatic potential (charge density), see Sachs … Woolf from ~early 2000s on the correct way to handle periodicity
 - Area per lipid: The value for double-tailed phospholipids is generally larger than single-chain hydrocarbons in systems like SAMs (60 $\AA^2$/molecule as opposed to 30-40 $\AA^2$/molecule).
 - Deuterium P$_2$ NMR order parameter: This parameter describes the alignment of lipid constituent bonds with the global membrane normal, and will vary along the length of the lipid tail group chains, and between liquid crystalline and gel phase lipids.
    - Area and NMR parameters are tightly coupled.  Area fluctuates on slow timescale (10s-100s of ns).
    - Computing error bars is tricky.  See order_parameters tool in LOOS for one approach
 - Lipid lateral diffusivity: This measure is used to characterize lipid mobility, to gain insight into collective lipid motion timescales, and also potentially to discriminate between liquid crystalline and gel phase lipids.
    - Be careful of artifacts due to lipid molecules partially or fully jumping across periodic boxes; account for periodicity by reimaging where appropriate.  Major box size dependence -- Yeh and Hummer
 - Thermodynamics/mechanics: area compressibility modulus: Theoretically (e.g. based on polymer brush theories), the K$_A$ is sometimes related to the oil-water interfacial tension.
 - Frame (box) area → contour area: this will allow for a deconvolution of undulations from area compressibility; to remove the effect of finite length scale undulations, results can also be extrapolated to a zero-sized membrane, where undulations no longer exist
 - Mechanics: fluctuation spectrum and bending modulus. For a tensionless membrane, the large wavelength/small wavevector behavior of the undulation spectrum (describing the height-height correlations of the membrane continuum shape) should follow an inverse fourth power relation in the wavevector, with a constant of proportionality that contains the bending modulus.
 - Lipid coordination number: this can be used in phase transitions and coexistence to distinguish between fluid and gel phases, which have markedly different coordination numbers. (There are several other structural, thermodynamic, and dynamic techniques for detecting and characterizing phase transitions and coexistence outlined above.)
 - Lateral radial distribution functions to characterize segregation of multicomponent bilayers (e.g. does cholesterol segregate with lipid X vs Y)
 - System size effects (cf. Waheed and Edholm)
    - Area per lipid
    - Area compressibility modulus
    - Bending modulus: as shown through renormalization group theory, the membrane softens at larger length scales, resulting in a lower bending modulus (cf. Pelitti and Leibler).
    - Diffusion coefficient (Yeh and Hummer)
    - Frank Brown had an excellent paper in ~2015 with Pastor and a few others

Diffusivity/PBC effects:
http://aip.scitation.org/doi/10.1063/1.4932980

Mechanical:
http://www.sciencedirect.com/science/article/pii/S0009308415300190

