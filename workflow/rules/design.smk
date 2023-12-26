

# List design attributes
rule list_design:
    output:
        touch("resources/tmp.log")
    shell:
        # f"echo wt={config['wt']}, mut={config['mut']}"
        "echo 123"