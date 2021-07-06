ARG BASEIMAGE=ontresearch/base-workflow-image:v0.1.0
FROM $BASEIMAGE
ARG ENVFILE=environment.yaml

COPY $ENVFILE $HOME/environment.yaml
RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \
    && micromamba install -n base --file $HOME/environment.yaml \
    # run medaka to download some models needed for variant calling
    && medaka tools download_models \
        --models r103_fast_variant_g507 r103_hac_variant_g507 r103_prom_variant_g3210 r103_sup_variant_g507 \
                 r941_min_fast_variant_g507 r941_min_hac_variant_g507 r941_min_sup_variant_g507 \
                 r941_prom_fast_variant_g507 r941_prom_hac_variant_g507 r941_prom_sup_variant_g507 \
                 r941_prom_variant_g303 r941_prom_variant_g322 r941_prom_variant_g360 \
    && micromamba clean --all --yes \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME

USER $WF_UID
WORKDIR $HOME
