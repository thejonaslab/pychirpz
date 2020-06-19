from fabric.api import local, env, run, put, cd, task
from fabric.contrib import project
import fabric as fab


env.roledefs['m'] = ['jonas@c65']

@task
def deploy():
    """
    Uses rsync. Note directory trailing-slash behavior is anti-intuitive
    """

    tgt_dir = "/data/jonas/pychirpz/"

    local('git ls-tree --full-tree --name-only -r HEAD |grep -v ipynb >  .git-files-list')
    local("echo '.git-files-list' >> .git-files-list")
    project.rsync_project("%s" % tgt_dir, local_dir="./",
                          exclude=['*.pickle', "*.ipynb"], # NOTE WILL NOT PUSH IPYTHON NOTEBOOKS

                          extra_opts='--files-from=.git-files-list')
    
    
    # project.rsync_project("%s/suncast/" % tgt_dir, 
    #                       local_dir="../setup.sh")
    # project.rsync_project("%s/data/info/" % tgt_dir, 
    #                       local_dir="../data/info/uploadqueries")
    # project.rsync_project("%s/data/info/" % tgt_dir, 
    #                       local_dir="../data/info/upload")
    
    # if "m" in env.effective_roles:
    
    #     project.rsync_project("%s/suncast/" % tgt_dir, local_dir="./",
                              
    #                           extra_opts="--include='*/' --include='*.jar' --exclude='*'")

    
    # if "m" in env.effective_roles:
    #     project.rsync_project("/data/jonas/suncast/suncast/notebooks.suncast/", local_dir="./notebooks.suncast",
    #                           extra_opts="--include '*.ipynb' --include '*.pdf' --include '*.png'  --include '*.mp4' --include='*/' --exclude='*' ", 
    #                           upload=False)
        

